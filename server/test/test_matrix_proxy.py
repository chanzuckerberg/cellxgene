import unittest
import numpy as np
from server.app.util.matrix_proxy import MatrixProxyView, MatrixProxy


class NdArrayProxyView(MatrixProxyView):
    """
    Fake test class for matrix proxy - wraps ndarray
    """

    @classmethod
    def __supports__(cls):
        return ('numpy.ndarray',)


class MatrixProxyViewTest(unittest.TestCase):

    def test_ismatrixproxy(self):
        n = np.zeros((2, 4))
        mp = MatrixProxy.create(n)

        self.assertIsNotNone(n)
        self.assertIsNotNone(mp)
        self.assertTrue(isinstance(mp, NdArrayProxyView))
        self.assertFalse(MatrixProxy.ismatrixproxy(n))
        self.assertTrue(MatrixProxy.ismatrixproxy(mp))

    def test_params(self):
        n = np.arange(15, dtype=np.float32).reshape((3, 5))
        mp = MatrixProxy.create(n)

        self.assertTrue(MatrixProxy.ismatrixproxy(mp))
        self.assertEqual(mp.ndim, 2)
        self.assertEqual(mp.shape, (3, 5))
        self.assertEqual(mp.dtype, np.float32)

        mpt = mp.T
        self.assertTrue(MatrixProxy.ismatrixproxy(mpt))
        self.assertEqual(mpt.ndim, 2)
        self.assertEqual(mpt.shape, (5, 3))
        self.assertEqual(mpt.dtype, np.float32)

    def test_toarray(self):
        n = np.arange(15, dtype=np.float32).reshape((3, 5))
        mp = MatrixProxy.create(n)
        self.assertTrue(
            np.all(mp.toarray() == [[0., 1., 2., 3., 4.], [5., 6., 7., 8., 9.],
                                    [10., 11., 12., 13., 14.]]))
        self.assertTrue(
            np.all(
                mp.T.toarray() == [[0., 5., 10.], [1., 6., 11.], [2., 7., 12.],
                                   [3., 8., 13.], [4., 9., 14.]]))

    def test_indexing(self):
        """
        [[ 0.,  1.,  2.,  3.,  4.],
         [ 5.,  6.,  7.,  8.,  9.],
         [10., 11., 12., 13., 14.]]
        """
        n = np.arange(15, dtype=np.float32).reshape((3, 5))
        mp = MatrixProxy.create(n)

        # should support int or slice for one or both dimensions

        # int, int

        self.assertEqual(mp[0, 0], 0)
        self.assertEqual(mp.T[0, 0], 0)
        self.assertEqual(mp[2, 4], 14)
        self.assertEqual(mp.T[4, 2], 14)
        self.assertEqual(mp[1, 3], mp.T[3, 1])

        # int, slice

        self.assertTrue(np.all(mp[0].toarray() == [0, 1, 2, 3, 4]))
        self.assertTrue(np.all(mp[1, 1:].toarray() == [6, 7, 8, 9]))
        self.assertTrue(np.all(mp[2, 1::-1].toarray() == [11, 10]))

        self.assertTrue(np.all(mp.T[0].toarray() == [0, 5, 10]))
        self.assertTrue(np.all(mp.T[1, 1:].toarray() == [6, 11]))
        self.assertTrue(np.all(mp.T[2, 1::-1].toarray() == [7, 2]))

        # slice, int

        self.assertTrue(np.all(mp[:, 0].toarray() == [0, 5, 10]))
        self.assertTrue(np.all(mp[1:, 1].toarray() == [6, 11]))
        self.assertTrue(np.all(mp[1::-1, 2].toarray() == [7, 2]))

        self.assertTrue(np.all(mp.T[:, 0].toarray() == [0, 1, 2, 3, 4]))
        self.assertTrue(np.all(mp.T[1:, 1].toarray() == [6, 7, 8, 9]))
        self.assertTrue(np.all(mp.T[1::-1, 2].toarray() == [11, 10]))

        # slice, slice

        self.assertTrue(np.all(mp[1:3, 2:4].toarray() == [[7, 8], [12, 13]]))
        self.assertTrue(
            np.all(mp[:3, :4].toarray() == [[0., 1., 2., 3.], [5., 6., 7., 8.],
                                            [10., 11., 12., 13.]]))
        self.assertTrue(
            np.all(mp[::-1, ::-1].toarray() ==
                   [[14, 13, 12, 11, 10], [9, 8, 7, 6, 5], [4, 3, 2, 1, 0]]))
        self.assertTrue(
            np.all(mp[::-2, ::-2].toarray() == [[14, 12, 10], [4, 2, 0]]))

        self.assertTrue(np.all(mp.T[2:4, 1:3].toarray() == [[7, 12], [8, 13]]))
        self.assertTrue(
            np.all(mp.T[:4, :3].toarray() == [[0, 5, 10], [1, 6, 11],
                                              [2, 7, 12], [3, 8, 13]]))
        self.assertTrue(
            np.all(mp.T[::-1, ::-1].toarray(
            ) == [[14, 9, 4], [13, 8, 3], [12, 7, 2], [11, 6, 1], [10, 5, 0]]))
        self.assertTrue(
            np.all(mp.T[::-2, ::-2].toarray() == [[14, 4], [12, 2], [10, 0]]))

    def test_repeated_indexing(self):
        """
        [[ 0.,  1.,  2.,  3.,  4.],
         [ 5.,  6.,  7.,  8.,  9.],
         [10., 11., 12., 13., 14.]]
        """
        n = np.arange(15, dtype=np.float32).reshape((3, 5))
        mp = MatrixProxy.create(n)

        self.assertEqual(mp[0][1], 1)
        self.assertEqual(mp.T[0][1], 5)

        self.assertTrue(np.all(mp[0::-1, ::-1][0, 2:4].toarray() == [2, 1]))
        self.assertTrue(np.all(mp[0::-1, 1:5:1][0, 1:3:1].toarray() == [2, 3]))
        self.assertTrue(np.all(mp[0::-1, 1:5:1][0, 2:0:-1].toarray() == [3, 2]))
        self.assertTrue(np.all(mp[0::-1, 5:1:-1][0, 1:3:1].toarray() == [3, 2]))
        self.assertTrue(np.all(mp[0::-1, 5:1:-1][0,
                                                 2:0:-1].toarray() == [2, 3]))

        self.assertTrue(np.all(mp.T[::-1, 0::-1][2:4, 0].toarray() == [2, 1]))
        self.assertTrue(np.all(mp.T[0::-1, 1:5:1][0, 1:3:1].toarray() == [10]))
        self.assertTrue(np.all(mp.T[0::-1, 1:5:1][0, 2:0:-1].toarray() == [10]))
        self.assertTrue(np.all(mp.T[0::-1, 5:1:-1][0, 1:3:1].toarray() == []))
        self.assertTrue(np.all(mp.T[0::-1, 5:1:-1][0, 2:0:-1].toarray() == []))

    def test_dimension_drop(self):
        """
        [[ 0.,  1.,  2.,  3.,  4.],
         [ 5.,  6.,  7.,  8.,  9.],
         [10., 11., 12., 13., 14.]]
        """
        n = np.arange(15, dtype=np.float32).reshape((3, 5))
        mp = MatrixProxy.create(n)

        # drop both dimensions, to a scalar
        self.assertEqual(mp[0, 0], 0)

        # drop 1 dimension, to an array
        self.assertTrue(np.all(mp[0, :].toarray() == [0, 1, 2, 3, 4]))
        self.assertTrue(np.all(mp[:, 0].toarray() == [0, 5, 10]))

        # with .T
        self.assertTrue(np.all(mp[0:2].T[-1:].toarray() == [[4, 9]]))
        self.assertTrue(np.all(mp[0:2].T[-1].toarray() == [4, 9]))

    def test_iter(self):
        """
        check that __iter__ is doing what we expect
        """
        n = np.arange(15, dtype=np.float32).reshape((3, 5))
        mp = MatrixProxy.create(n)

        rows = [r for r in mp]
        self.assertEqual(len(rows), 3)
        self.assertTrue(np.all(rows[0].toarray() == [0, 1, 2, 3, 4]))
        for i, r in enumerate(rows):
            self.assertTrue(np.all(mp[i].toarray() == r.toarray()))

        cols = [c for c in mp.T]
        self.assertEqual(len(cols), 5)
        self.assertTrue(np.all(cols[0].toarray() == [0, 5, 10]))
        for i, c in enumerate(cols):
            self.assertTrue(np.all(mp.T[i].toarray() == c.toarray()))

        e = [e for e in mp[0]]
        self.assertEqual(len(e), 5)
        self.assertEqual(e, [0, 1, 2, 3, 4])
