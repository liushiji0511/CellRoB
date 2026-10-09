import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import os
import math

class CellRoB:
    def __init__(self, **kwargs):
        self.filelist_dir = kwargs.get('filelist_dir', [])
        self.binsize = kwargs.get('binsize', 40000)
        self.result_dir = kwargs.get('result_dir', 'cellrob_result')
        self.chroms = kwargs.get('chroms', 'chr6')
        self.chromsizes_dir = kwargs.get('chromsizes_dir')
        self.impute_n = kwargs.get('impute_n', 15)
        self.window_len = kwargs.get('window', 5)
        self.Knn_path = kwargs.get('Knn_path')
        self.worker_n = kwargs.get('worker_n', 1)

        self.filelist = pd.read_csv(self.filelist_dir, header=None).squeeze().tolist()
        self.Knn = np.load(self.Knn_path)
        self.chromsizes = pd.read_csv(self.chromsizes_dir, header=None, sep='\t', names=['chrom', 'size'])

        self.breakpoint = 40000 // self.binsize
        self.MAX_DIS = 500 * 40000 // self.binsize

        print(f'num of file: {len(self.filelist)}')
        print(f'result dir: {self.result_dir}')
        print(f'binsize: {self.binsize}; window size: {self.window_len}; impute_n: {self.impute_n}')

    def impute(self, filepath):
        contact = pd.read_csv(filepath, sep='\t', header=0)
        contact['count'] = 1
        cell_index = self.filelist.index(filepath)
        knn = self.Knn[cell_index].copy()
        knn[cell_index] = 0
        knn[~np.isfinite(knn)] = 0
        for _ in range(self.impute_n):
            maxindex = np.argmax(knn)
            similarity = knn[maxindex]
            if similarity <= 0.05:
                break
            knn[maxindex] = 0
            simcontact = pd.read_csv(self.filelist[maxindex], sep='\t', header=0)
            simcontact['count'] = similarity
            contact = pd.concat((contact, simcontact), axis=0)
        return contact

    def matrix_to_insulationScore(self, mat):
        _len = mat.shape[0]
        windowlen = self.window_len
        IS = np.zeros(mat.shape[0], dtype=float)
        for i in range(windowlen, _len - windowlen):
            c = mat[i - windowlen:i, i + 1:i + windowlen].sum()
            a = 0
            k = i - windowlen
            while k < i:
                a += mat[k, k + 1:i].sum()
                k += 1
            b = 0
            k = i + 1
            while k < i + windowlen:
                b += mat[k, k + 1:i + windowlen + 1].sum()
                k += 1
            if a + b - c <= 0:
                IS[i] = 0
            else:
                IS[i] = (a + b - c) / (a + b + c)
        return IS

    def detect_can_boundaries(self, vector):
        MAX = 1e6
        max_dis = self.MAX_DIS
        hdmd_vector = np.empty(vector.shape[0])
        vector = vector + (np.random.rand(vector.shape[0]) / MAX)
        for loc in range(vector.shape[0]):
            element = vector[loc]
            if element != 1:
                upper_limit = max(-1, loc - max_dis - 1)
                flag1 = 0
                for left in range(loc - 1, upper_limit, -1):
                    if vector[left] > element:
                        flag1 = 1
                        break
                flag2 = 0
                lower_limit = min(vector.shape[0], loc + max_dis + 1)

                for right in range(loc + 1, lower_limit):
                    if vector[right] > element:
                        flag2 = 1
                        break
                if (flag1 == 0) and (flag2 != 0):
                    hdmdp = right - loc
                elif (flag1 != 0) and (flag2 == 0):
                    hdmdp = loc - left
                elif (flag1 != 0) and (flag2 != 0):
                    hdmdp = min(loc - left, right - loc)
                else:
                    hdmdp = max_dis
                hdmd_vector[loc] = hdmdp
            else:
                hdmd_vector[loc] = np.nan
        return hdmd_vector

    def cluster_result(self, boundaries_record):
        boundaries_record['den*hdmd'] = boundaries_record['den'] * boundaries_record['hdmd']

        boundaries_record['rank'] = boundaries_record['den*hdmd'].rank(ascending=False, method='first')
        temp_rank = boundaries_record['rank'] / np.nanmax(boundaries_record['rank'])

        temp_eta = boundaries_record['den*hdmd'] / np.nanmax(boundaries_record['den*hdmd'])

        boundaries_record['transformed_den*hdmd'] = (temp_eta + temp_rank) / np.sqrt(2)  # 可以对boundaries进行更好的区分
        breakpoint = boundaries_record.iloc[boundaries_record['transformed_den*hdmd'].idxmin()]['den*hdmd']

        boundaries_record['cluster'] = -1

        boundaries_record['cluster'] = boundaries_record['den*hdmd'] / breakpoint

        return boundaries_record

    def mat_to_bd(self, mat):

        np.fill_diagonal(mat, 0)
        mat = mat + mat.T
        bin_quality = np.ones(len(mat), dtype=np.int64)
        bin_quality[np.sum(mat, axis=0) < 1] = 0

        if np.sum(bin_quality) < 100:
            return np.zeros(len(mat)), np.zeros(len(mat))

        nice_mat = mat[:, bin_quality == 1][bin_quality == 1, :]

        IS = self.matrix_to_insulationScore(nice_mat)
        hdmd_IS = self.detect_can_boundaries(IS)
        boundaries_record = pd.DataFrame()
        boundaries_record['den'] = IS
        boundaries_record['hdmd'] = hdmd_IS - 1

        boundaries_record = self.cluster_result(boundaries_record)
        bd_rank = boundaries_record['cluster'].values


        full_IS = np.zeros(len(bin_quality))
        full_bd = np.zeros(len(bin_quality))
        full_IS[bin_quality == 1] = IS
        full_bd[bin_quality == 1] = bd_rank

        return full_bd.flatten(), full_IS.flatten()

    def run(self):
        def process_file(filepath):
            filename = filepath.split('/')[-1]
            print(f'processing: {filename}')

            if os.path.exists(self.result_dir + f'/{self.chroms[-1]}/{filename}.boundary.npy'):
                return True

            contact = self.impute(filepath)
            print('impute done')
            contact[['pos1', 'pos2']] = contact[['pos1', 'pos2']] // self.binsize
            contact = contact.groupby(['chrom1', 'pos1', 'chrom2', 'pos2'])['count'].sum().reset_index()
            print('binned')

            for chrom in self.chroms:

                chrom_contact = contact[(contact['chrom1'] == chrom) & (contact['chrom2'] == chrom)].copy()
                chrom_contact = chrom_contact[['pos1', 'pos2', 'count']].values
                _len = self.chromsizes['size'][self.chromsizes['chrom'] == chrom].values[0] // self.binsize + 1
                mat = np.zeros((_len, _len))

                for row in chrom_contact:
                    mat[int(row[0]), int(row[1])] += row[2]
                    mat[int(row[1]), int(row[0])] += row[2]  # 确保对称性

                bd_rank, insulationScore = self.mat_to_bd(mat)

                os.makedirs(os.path.join(self.result_dir, chrom), exist_ok=True)

                np.save(os.path.join(self.result_dir, chrom, f'{filename}.boundary.npy'), bd_rank)
                np.save(os.path.join(self.result_dir, chrom, f'{filename}.insulation.npy'), insulationScore)


        with ThreadPoolExecutor(max_workers=self.worker_n) as executor:
            executor.map(process_file, self.filelist)

        print('hi')

    def concatenate_boundary_insulationScore(self):
        filenames = [filename.split('/')[-1] for filename in self.filelist]
        for chrom in self.chroms:
            print(f'concatenate: {chrom}')
            bds_re = []
            ins_re = []
            for filename in filenames:
                bd = np.load(os.path.join(self.result_dir, chrom, f'{filename}.boundary.npy'))
                ins = np.load(os.path.join(self.result_dir, chrom, f'{filename}.insulation.npy'))
                bds_re.append(bd.reshape((1, len(bd))))
                ins_re.append(ins.reshape((1, len(ins))))
            bds = np.concatenate(bds_re)
            ins = np.concatenate(ins_re)
            os.makedirs(os.path.join(self.result_dir, 'concatenate'), exist_ok=True)
            np.save(os.path.join(self.result_dir, 'concatenate', f'{chrom}.boundary.npy'), bds)
            np.save(os.path.join(self.result_dir, 'concatenate', f'{chrom}.insulation.npy'), ins)

    def align(self, count_threshold=0):
        def subdivide(signal):
            from scipy.signal import find_peaks

            boundaries = find_peaks(-signal, height=None, threshold=None,
                                    prominence=None, width=None, rel_height=0.5)[0]
            _boundaries = []
            for index in range(len(boundaries)):
                if signal[boundaries[index] - 1] > signal[boundaries[index] + 1]:
                    _boundaries.append(boundaries[index] + 0.5)
                else:
                    _boundaries.append(boundaries[index] - 0.5)

            _boundaries = [0] + [i - 1 for i in _boundaries] + [len(signal) - 3]
            return _boundaries

        def get_merge_anchor(bds, IS):
            _len = bds.shape[1]
            loss_record = np.zeros(_len)

            for ship in range(_len):
                loss = 0
                for detector in range(_len):
                    loss += np.sum(IS[bds[:, detector] == 1, detector] - IS[bds[:, detector] == 1, ship])
                loss_record[ship] = loss

            anchor = np.argmin(loss_record)

            return anchor

        for chrom in self.chroms:
            bds = np.load(os.path.join(self.result_dir, 'concatenate', f'{chrom}.boundary.npy'))
            IS = np.load(os.path.join(self.result_dir, 'concatenate', f'{chrom}.insulation.npy'))
            IS_row_mean = np.mean(IS, axis=1, keepdims=True)
            bds[IS <= IS_row_mean] = 0
            bds = np.where(bds >= self.breakpoint, 1, 0)
            bv = np.sum(bds, axis=0)

            kernel_size = 40000 // self.binsize
            kernel = np.ones(kernel_size * 2 + 1)
            bv = np.convolve(bv, kernel, mode='same')

            signal = bv
            signal = np.concatenate(([signal[0] * 0.5], signal, [signal[-1] * 0.5]))
            boundaries = subdivide(signal)

            _regions = np.array([[0, 0]])
            for i in range(len(boundaries) - 1):
                _regions = np.concatenate(
                    (_regions, np.array([[math.ceil(boundaries[i]), math.floor(boundaries[i + 1])]])))
            regions = _regions[1:,:]

            record = np.zeros(bds.shape[1])
            Rob_bds = np.zeros_like(bds, dtype=np.float64)
            for region in regions:
                start = region[0]
                end = region[1]
                bds_slice = bds[:, start:end + 1]
                IS_slice = IS[:, start:end + 1]
                cell_with_bd = np.sum(bds_slice, axis=1) != 0
                if bds_slice.sum() > count_threshold:
                    if start - end != 0:
                        anchor = get_merge_anchor(bds_slice, IS_slice)
                        cen = start + anchor
                        Rob_bds[cell_with_bd, cen] = np.max(IS[cell_with_bd, start: end+1], axis=1)
                    else:
                        Rob_bds[cell_with_bd, start] = IS[cell_with_bd, start]
                    record[start:end + 1] += 1

            os.makedirs(os.path.join(self.result_dir, 'align'), exist_ok=True)
            np.save(os.path.join(self.result_dir, 'align', f'{chrom}.boundary.npy'), Rob_bds)
