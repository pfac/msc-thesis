typedef unsigned long ulong;

template<typename T>
void sqrtm (const Mat<T>& t, Mat<T>& u) {
	const ulong n = t.n_rows;
	u = zeros<Mat<T> >(n,n);
	for (ulong d = 0; d < n; ++d) {
		const ulong m = n - d;
#pragma omp parallel for
		for (ulong e = 0; e < m; ++e) {
			if (d == 0)
				u(e,e) = std::sqrt( t(e,e) );
			else {
				// find the index
				const ulong i = e;
				const ulong j = e + d;
				// solve dependencies
				const ulong z[2] = { i + 1, j - 1 };
				T s = T(0);
				if (z[1] > z[0]) {
					const span k(z[0], z[1]);
					const Mat<T> tmp = r(i,k) * r(k,j);
					s = tmp(0,0);
				}
				r(i,j) = (t(i,j) - s) / (r(i,i) + r(j,j));
			}
		}
	}
}
