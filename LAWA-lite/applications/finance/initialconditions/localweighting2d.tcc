namespace lawa {


template <typename T, typename Basis>
T
LocalWeighting2D<T,Basis>::left_x1;

template <typename T, typename Basis>
T
LocalWeighting2D<T,Basis>::right_x1;

template <typename T, typename Basis>
T
LocalWeighting2D<T,Basis>::left_x2;

template <typename T, typename Basis>
T
LocalWeighting2D<T,Basis>::right_x2;

template <typename T, typename Basis>
T
LocalWeighting2D<T,Basis>::RightmLeft_x1;

template <typename T, typename Basis>
T
LocalWeighting2D<T,Basis>::RightmLeft_x2;

template <typename T, typename Basis>
const Basis *
LocalWeighting2D<T,Basis>::basis;

template <typename T, typename Basis>
int
LocalWeighting2D<T,Basis>::weight_type;

template <typename T, typename Basis>
void
LocalWeighting2D<T,Basis>::setDomain(T _left_x1, T _right_x1, T _left_x2, T _right_x2)
{
    left_x1  = _left_x1;
    right_x1 = _right_x1;
    RightmLeft_x1 = right_x1 - left_x1;

    left_x2  = _left_x2;
    right_x2 = _right_x2;
    RightmLeft_x2 = right_x2 - left_x2;
}

template <typename T, typename Basis>
void
LocalWeighting2D<T,Basis>::setBasis(const Basis *_basis)
{
    basis = _basis;
}

template <typename T, typename Basis>
void
LocalWeighting2D<T,Basis>::setWeightType(int _weight_type)
{
    weight_type = _weight_type;
}

template <typename T, typename Basis>
T
LocalWeighting2D<T,Basis>::weight(const Index2D &index)
{
    if (weight_type == 0) return 1.;

    int  j_x=index.index1.j,    j_y=index.index2.j;
    long k_x=index.index1.k,    k_y=index.index2.k;
    XType xtype_x=index.index1.xtype, xtype_y=index.index2.xtype;


    Support<T> supp_x = (*basis).first.generator(xtype_x).support(j_x,k_x);
    T center_x  = 0.5*(supp_x.l1 + supp_x.l2);
    Support<T> supp_y = (*basis).second.generator(xtype_y).support(j_y,k_y);
    T center_y  = 0.5*(supp_y.l1 + supp_y.l2);
    if ((Basis::FirstBasisType::Domain!=R) && (Basis::SecondBasisType::Domain!=R)) {
        center_x *= RightmLeft_x1; center_x += left_x1;
        center_y *= RightmLeft_x1; center_y += left_x1;
    }
    //std::cerr << index << " " << supp_x << " " << supp_y << ":  (" << center_x << ", " << center_y << ")" << std::endl;

    if (weight_type == 1) {
        return std::exp(-3.*(center_x*center_x + center_y*center_y));
    }
    else if (weight_type == 2) {    // only for realline!!
        T x1InnerMax = max(fabs(left_x1),fabs(right_x1))-1.;
        T x2InnerMax = max(fabs(left_x2),fabs(right_x2))-1.;
        T dist_center_x = fabs(fabs(center_x)-x1InnerMax);
        T dist_center_y = fabs(fabs(center_y)-x2InnerMax);
        if (fabs(center_x)<=x1InnerMax && fabs(center_y)<=x2InnerMax) {
            return 1.;
        }
        else {
            std::cerr << index << " -> " << center_x << " " << dist_center_x << " :: "
                      << center_y << " " << dist_center_y << std::endl;
            return std::exp(-3.*(dist_center_x*dist_center_x + dist_center_y*dist_center_y));
        }
    }
}

}   // namespace lawa
