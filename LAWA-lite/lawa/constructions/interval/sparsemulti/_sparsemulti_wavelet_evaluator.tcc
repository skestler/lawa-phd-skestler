namespace lawa {
//--- cubic evaluators ---------------------------------------------------------

template <typename T>
T
_sparsemulti_cubic_wavelet_inner_evaluator0(T x, unsigned short deriv)
{
    T value = 0.0;
    if (0.<=x && x<=2.) {
        value =
        -(2./15.)*_cubic_sparsemulti_scaling_inner_evaluator0(2*x-1,deriv)
        +(4./15.)*_cubic_sparsemulti_scaling_inner_evaluator0(2*x-2,deriv)
        -(2./15.)*_cubic_sparsemulti_scaling_inner_evaluator0(2*x-3,deriv)

        -         _cubic_sparsemulti_scaling_inner_evaluator1(2*x-1,deriv)
        +         _cubic_sparsemulti_scaling_inner_evaluator1(2*x-3,deriv);

        value *= pow2i<T>(deriv);
    }
    return value;
}

template <typename T>
T
_sparsemulti_cubic_wavelet_inner_evaluator1(T x, unsigned short deriv)
{
    T value = 0.0;

    if (0.<=x && x<=2.) {
        value =
         (7./39.)* _cubic_sparsemulti_scaling_inner_evaluator0(2*x-1,deriv)

        -(7./39.)* _cubic_sparsemulti_scaling_inner_evaluator0(2*x-3,deriv)

        +          _cubic_sparsemulti_scaling_inner_evaluator1(2*x-1,deriv)
        +(44./13.)*_cubic_sparsemulti_scaling_inner_evaluator1(2*x-2,deriv)
        +          _cubic_sparsemulti_scaling_inner_evaluator1(2*x-3,deriv);

        value *= pow2i<T>(deriv);
    }
    return value;
}

template <typename T>
T
_sparsemulti_cubic_wavelet_inner_evaluator2(T x, unsigned short deriv)
{
    T value = 0.0;
    
    if (-2.<=x && x<=2.) {
        value =
        -(4595./13728.)  *_cubic_sparsemulti_scaling_inner_evaluator0(2*x+3,deriv)
        +(7./65.)        *_cubic_sparsemulti_scaling_inner_evaluator0(2*x+2,deriv)
        -(18737./68640.) *_cubic_sparsemulti_scaling_inner_evaluator0(2*x+1,deriv)
        +                 _cubic_sparsemulti_scaling_inner_evaluator0(2*x  ,deriv)
        -(18737./68640.) *_cubic_sparsemulti_scaling_inner_evaluator0(2*x-1,deriv)
        +(7./65.)        *_cubic_sparsemulti_scaling_inner_evaluator0(2*x-2,deriv)
        -(4595./13728.)  *_cubic_sparsemulti_scaling_inner_evaluator0(2*x-3,deriv)

        -(68741./22880.) *_cubic_sparsemulti_scaling_inner_evaluator1(2*x+3,deriv)
        -(69./40.)       *_cubic_sparsemulti_scaling_inner_evaluator1(2*x+2,deriv)
        -(204701./22880.)*_cubic_sparsemulti_scaling_inner_evaluator1(2*x+1,deriv)

        +(204701./22880.)*_cubic_sparsemulti_scaling_inner_evaluator1(2*x-1,deriv)
        +(69./40.)       *_cubic_sparsemulti_scaling_inner_evaluator1(2*x-2,deriv)
        +(68741./22880.) *_cubic_sparsemulti_scaling_inner_evaluator1(2*x-3,deriv);

        value *= pow2i<T>(deriv);
    }
    return value;
}

template <typename T>
T
_sparsemulti_cubic_wavelet_inner_evaluator3(T x, unsigned short deriv)
{
    T value = 0.0;

    if (-2.<=x && x<=2.) {
        value =
        +(417./22880.)   *_cubic_sparsemulti_scaling_inner_evaluator0(2*x+3,deriv)
        -(7./2340.)      *_cubic_sparsemulti_scaling_inner_evaluator0(2*x+2,deriv)
        +(5443./205920.) *_cubic_sparsemulti_scaling_inner_evaluator0(2*x+1,deriv)

        -(5443./205920.) *_cubic_sparsemulti_scaling_inner_evaluator0(2*x-1,deriv)
        +(7./2340.)      *_cubic_sparsemulti_scaling_inner_evaluator0(2*x-2,deriv)
        -(417./22880.)   *_cubic_sparsemulti_scaling_inner_evaluator0(2*x-3,deriv)

        +(723./4576.)    *_cubic_sparsemulti_scaling_inner_evaluator1(2*x+3,deriv)
        +(1./8.)         *_cubic_sparsemulti_scaling_inner_evaluator1(2*x+2,deriv)
        +(8153./13728.)  *_cubic_sparsemulti_scaling_inner_evaluator1(2*x+1,deriv)
        +(1./2.)         *_cubic_sparsemulti_scaling_inner_evaluator1(2*x  ,deriv)
        +(8153./13728.)  *_cubic_sparsemulti_scaling_inner_evaluator1(2*x-1,deriv)
        +(1./8.)         *_cubic_sparsemulti_scaling_inner_evaluator1(2*x-2,deriv)
        +(723./4576.)    *_cubic_sparsemulti_scaling_inner_evaluator1(2*x-3,deriv);

        value *= pow2i<T>(deriv);
    }
    return value;
}

template <typename T>
T
_sparsemulti_cubic_wavelet_left_evaluator0(T x, unsigned short deriv)
{
    return _sparsemulti_cubic_wavelet_inner_evaluator0(x, deriv);
}

template <typename T>
T
_sparsemulti_cubic_wavelet_left_evaluator1(T x, unsigned short deriv)
{
    return _sparsemulti_cubic_wavelet_inner_evaluator1(x, deriv);
}

template <typename T>
T
_sparsemulti_cubic_wavelet_left_evaluator2(T x, unsigned short deriv)
{
    T value = 0.0;

    if (0<=x && x<=2.) {
        value =
        +(417./22880.)   *_cubic_sparsemulti_scaling_inner_evaluator0(2*x+3,deriv)
        -(7./2340.)      *_cubic_sparsemulti_scaling_inner_evaluator0(2*x+2,deriv)
        +(5443./205920.) *_cubic_sparsemulti_scaling_inner_evaluator0(2*x+1,deriv)

        -(5443./205920.) *_cubic_sparsemulti_scaling_inner_evaluator0(2*x-1,deriv)
        +(7./2340.)      *_cubic_sparsemulti_scaling_inner_evaluator0(2*x-2,deriv)
        -(417./22880.)   *_cubic_sparsemulti_scaling_inner_evaluator0(2*x-3,deriv)

        +(723./4576.)    *_cubic_sparsemulti_scaling_inner_evaluator1(2*x+3,deriv)
        +(1./8.)         *_cubic_sparsemulti_scaling_inner_evaluator1(2*x+2,deriv)
        +(8153./13728.)  *_cubic_sparsemulti_scaling_inner_evaluator1(2*x+1,deriv)
        +(1./2.)         *_cubic_sparsemulti_scaling_inner_evaluator1(2*x  ,deriv)
        +(8153./13728.)  *_cubic_sparsemulti_scaling_inner_evaluator1(2*x-1,deriv)
        +(1./8.)         *_cubic_sparsemulti_scaling_inner_evaluator1(2*x-2,deriv)
        +(723./4576.)    *_cubic_sparsemulti_scaling_inner_evaluator1(2*x-3,deriv);

        value *= pow2i<T>(deriv);
    }
    return value;
}
    
template <typename T>
T
_sparsemulti_cubic_wavelet_right_evaluator0(T x, unsigned short deriv)
{
    T value = 0.0;

    if (-2.<=x && x<=0.) {
        value =
        +(417./22880.)   *_cubic_sparsemulti_scaling_inner_evaluator0(2*x+3,deriv)
        -(7./2340.)      *_cubic_sparsemulti_scaling_inner_evaluator0(2*x+2,deriv)
        +(5443./205920.) *_cubic_sparsemulti_scaling_inner_evaluator0(2*x+1,deriv)

        -(5443./205920.) *_cubic_sparsemulti_scaling_inner_evaluator0(2*x-1,deriv)
        +(7./2340.)      *_cubic_sparsemulti_scaling_inner_evaluator0(2*x-2,deriv)
        -(417./22880.)   *_cubic_sparsemulti_scaling_inner_evaluator0(2*x-3,deriv)

        +(723./4576.)    *_cubic_sparsemulti_scaling_inner_evaluator1(2*x+3,deriv)
        +(1./8.)         *_cubic_sparsemulti_scaling_inner_evaluator1(2*x+2,deriv)
        +(8153./13728.)  *_cubic_sparsemulti_scaling_inner_evaluator1(2*x+1,deriv)
        +(1./2.)         *_cubic_sparsemulti_scaling_inner_evaluator1(2*x  ,deriv)
        +(8153./13728.)  *_cubic_sparsemulti_scaling_inner_evaluator1(2*x-1,deriv)
        +(1./8.)         *_cubic_sparsemulti_scaling_inner_evaluator1(2*x-2,deriv)
        +(723./4576.)    *_cubic_sparsemulti_scaling_inner_evaluator1(2*x-3,deriv);

        value *= pow2i<T>(deriv);
    }
    return value;
}


} // namespace lawa
