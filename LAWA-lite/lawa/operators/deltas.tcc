namespace lawa {

template <typename T, typename Basis>
GeMatrix<FullStorage<T,ColMajor> >
computeDeltas(const Basis &basis, int j, long k, XType e)
{
    GeMatrix<FullStorage<T,ColMajor> > ret;
    const typename Basis::BasisFunctionType &varphi=basis.generator(e); //either B-Spline or Wavelet
    int d = basis.d;

    ret.engine().resize(varphi.singularSupport(j,k).length(),2);
    ret(_,1) = varphi.singularSupport(j,k);
    Support<T> supp = varphi.support(j,k);
    T step = pow2i<T>(-(j+5)); // 1.0/(1<<(j+1));
    for (int i = 1; i<=ret.numRows(); ++i) {
        ret(i,2) = ((ret(i,1)==supp.l2) ? 0.0 : varphi(std::min(ret(i,1)+step, supp.l2),j,k,d-1))
                 - ((ret(i,1)==supp.l1) ? 0.0 : varphi(std::max(ret(i,1)-step, supp.l1),j,k,d-1));
    }
    return ret;
}


}    //namespace lawa
