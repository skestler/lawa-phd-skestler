namespace lawa {

template <typename T>
T
TruncatedBasketPutOption2D<T>::left_x1;

template <typename T>
T
TruncatedBasketPutOption2D<T>::right_x1;

template <typename T>
T
TruncatedBasketPutOption2D<T>::left_x2;

template <typename T>
T
TruncatedBasketPutOption2D<T>::right_x2;

template <typename T>
flens::DenseVector<Array<T> >
TruncatedBasketPutOption2D<T>::singPts_x1;

template <typename T>
flens::DenseVector<Array<T> >
TruncatedBasketPutOption2D<T>::singPts_x2;


template <typename T>
Option2D<T,BasketPut>
TruncatedBasketPutOption2D<T>::basketputoption;



template <typename T>
T
TruncatedBasketPutOption2D<T>::u11;

template <typename T>
T
TruncatedBasketPutOption2D<T>::u21;

template <typename T>
T
TruncatedBasketPutOption2D<T>::u12;

template <typename T>
T
TruncatedBasketPutOption2D<T>::u22;


template <typename T>
int
TruncatedBasketPutOption2D<T>::type;

template <typename T>
T
TruncatedBasketPutOption2D<T>::truncWidth;

template <typename T>
T
TruncatedBasketPutOption2D<T>::damping_c;


template <typename T>
DenseVector<Array<T> >
TruncatedBasketPutOption2D<T>::critical_line_x1;

template <typename T>
bool
TruncatedBasketPutOption2D<T>::critical_above_x1;

template <typename T>
DenseVector<Array<T> >
TruncatedBasketPutOption2D<T>::critical_line_x2;

template <typename T>
bool
TruncatedBasketPutOption2D<T>::critical_above_x2;

template <typename T>
void
TruncatedBasketPutOption2D<T>::setOption(const Option2D<T,BasketPut> &_basketputoption)
{
    basketputoption.optionparameters = _basketputoption.optionparameters;
    u11 = 1.;
    u21 = 0.;
    u12 = 0.;
    u22 = 1.;
}

template <typename T>
void
TruncatedBasketPutOption2D<T>::setTransformation(T _u11, T _u21, T _u12, T _u22)
{
    u11 = _u11;
    u21 = _u21;
    u12 = _u12;
    u22 = _u22;
}

template <typename T>
void
TruncatedBasketPutOption2D<T>::setTruncation(T _left_x1, T _right_x1, T _left_x2, T _right_x2,
                                             int _type, T _truncWidth, T _damping_c)
{
    left_x1      = _left_x1;
    right_x1     = _right_x1;
    left_x2      = _left_x2;
    right_x2     = _right_x2;
    type         = _type;
    truncWidth   = _truncWidth;
    damping_c    = _damping_c;

//    singPts_x1.engine().resize(10);
//    singPts_x1(1) = left_x1 + truncWidth, -0.04, -0.03, -0.02, -0.01, 0., 0.01, 0.02, 0.03, 0.04;

//    singPts_x2.engine().resize(10);
//    singPts_x2(1) = left_x2 + truncWidth, -0.04, -0.03, -0.02, -0.01, 0., 0.01, 0.02, 0.03, 0.04;
    singPts_x1.engine().resize(1);
    singPts_x1(1) = left_x1 + truncWidth;

    singPts_x2.engine().resize(1);
    singPts_x2(1) = left_x2 + truncWidth;
}

template <typename T>
void
TruncatedBasketPutOption2D<T>::setCriticalLine_x1(T _critical_line_x1, bool _critical_above_line_x1)
{
    critical_line_x1.engine().resize(1);
    critical_line_x1(1) = _critical_line_x1;
    critical_above_x1   = _critical_above_line_x1;
}

template <typename T>
void
TruncatedBasketPutOption2D<T>::setCriticalLine_x2(T _critical_line_x2, bool _critical_above_line_x2)
{
    critical_line_x2.engine().resize(1);
    critical_line_x2(1) = _critical_line_x2;
    critical_above_x2   = _critical_above_line_x2;
}

template <typename T>
bool
TruncatedBasketPutOption2D<T>::isCritical(T a1, T b1, T a2, T b2)
{
    if (critical_line_x1.length()>0) {
        if (critical_above_x1) {
            if (!(a1>critical_line_x1(1))) return true;
        }
        else {
            if (!(b1<critical_line_x1(1))) return true;
        }
    }
    if (critical_line_x2.length()>0) {
        if (critical_above_x2) {
            if (!(a2>critical_line_x2(1))) return true;
        }
        else {
            if (!(b2<critical_line_x2(1))) return true;
        }
    }

    return false;
}


template <typename T>
T
TruncatedBasketPutOption2D<T>::payoff(T x1, T x2)
{
    if (x1>=left_x1+truncWidth && x2>=left_x2+truncWidth) {
        return basketputoption.payoff_log(u11*x1 + u21*x2, u12*x1 + u22*x2);
    }
    else if (x1<left_x1+truncWidth && x2>=left_x2+truncWidth) {
        T dist1 = (left_x1+truncWidth) - x1;
        return exp(-damping_c*dist1*dist1)*basketputoption.payoff_log(u11*x1 + u21*x2, u12*x1 + u22*x2);
    }
    else if (x1>=left_x1+truncWidth && x2<left_x2+truncWidth) {
        T dist2 = (left_x2+truncWidth) - x2;
        return exp(-damping_c*dist2*dist2)*basketputoption.payoff_log(u11*x1 + u21*x2, u12*x1 + u22*x2);
    }
    else {
        T dist1 = (left_x1+truncWidth) - x1;
        T dist2 = (left_x2+truncWidth) - x2;
        return exp(-damping_c*dist1*dist1)*exp(-damping_c*dist2*dist2)*basketputoption.payoff_log(u11*x1 + u21*x2, u12*x1 + u22*x2);
    }

    return basketputoption.payoff_log(u11*x1 + u21*x2, u12*x1 + u22*x2);
}

}   // namespace lawa
