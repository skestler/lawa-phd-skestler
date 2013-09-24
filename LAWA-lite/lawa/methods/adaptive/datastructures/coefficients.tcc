/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Kristina Steih, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

namespace lawa {

/* **********************************************************************************************
 * Lexicographically sorted coefficient vector
 * ********************************************************************************************** */

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>::Coefficients(void)
{
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>&
Coefficients<Lexicographical,T,Index>::operator=(const Coefficients<Lexicographical,T,Index> &_coeff)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    typedef typename Coefficients<Lexicographical,T,Index>::value_type val_type;
    this->erase(Coefficients<Lexicographical,T,Index>::begin(), Coefficients<Lexicographical,T,Index>::end());

    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            this->insert(val_type((*lambda).first, (*lambda).second));
        }
    }

    return *this;
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>&
Coefficients<Lexicographical,T,Index>::operator=(const Coefficients<AbsoluteValue,T,Index> &_coeff)
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_it;
    typedef typename Coefficients<Lexicographical,T,Index>::value_type val_type;
    this->erase(Coefficients<Lexicographical,T,Index>::begin(), Coefficients<Lexicographical,T,Index>::end());

    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            this->insert(val_type((*lambda).second, (*lambda).first));
        }
    }

    return *this;
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>
Coefficients<Lexicographical,T,Index>::operator-(const Coefficients<Lexicographical,T,Index> &_coeff) const
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    Coefficients<Lexicographical,T,Index> ret = *this;
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            ret.operator[]((*lambda).first) -= (*lambda).second;
        }
    }
    return ret;
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index> &
Coefficients<Lexicographical,T,Index>::operator-=(const Coefficients<Lexicographical,T,Index> &_coeff)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    //Coefficients<Lexicographical,T,Index> ret = *this;
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            (*this).operator[]((*lambda).first) -= (*lambda).second;
        }
    }
    return (*this);
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index> &
Coefficients<Lexicographical,T,Index>::operator+=(const Coefficients<Lexicographical,T,Index> &_coeff)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    //Coefficients<Lexicographical,T,Index> ret = *this;
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            (*this).operator[]((*lambda).first) += (*lambda).second;
        }
    }
    return (*this);
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index> &
Coefficients<Lexicographical,T,Index>::operator*=(const T factor)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    //Coefficients<Lexicographical,T,Index> ret = *this;
    if ((*this).size() > 0) {
        for (const_it lambda = (*this).begin(); lambda != (*this).end(); ++lambda) {
            (*this).operator[]((*lambda).first) *= factor;
        }
    }
    return (*this);
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>
Coefficients<Lexicographical,T,Index>::operator+(const Coefficients<Lexicographical,T,Index> &_coeff) const
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    Coefficients<Lexicographical,T,Index> ret = *this;
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            ret.operator[]((*lambda).first) += (*lambda).second;
        }
    }
    return ret;
}

template <typename T, typename Index>
T
Coefficients<Lexicographical,T,Index>::operator*(const Coefficients<Lexicographical,T,Index> &_coeff2) const
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    //Coefficients<Lexicographical,T,Index> _coeff1 = *this;
    long double ret = 0.L;
    if (_coeff2.size() > 0) {
        for (const_it lambda = _coeff2.begin(); lambda != _coeff2.end(); ++lambda) {
            const_it mu = (*this).find(((*lambda).first));
            if (mu!=(*this).end()) {
                ret += (long double)( (*mu).second * (*lambda).second );
            }
            //ret += (long double)((*this).operator[]((*lambda).first) * (*lambda).second);
        }
    }
    return (T)ret;
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>
Coefficients<Lexicographical,T,Index>::operator*(const T factor) const
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    Coefficients<Lexicographical,T,Index> ret = *this;
	for (const_it lambda = (*this).begin(); lambda != (*this).end(); ++lambda) {
        ret.operator[]((*lambda).first) *= factor;
    }
    
    return ret;
}

template <typename T, typename Index>
void
Coefficients<Lexicographical,T,Index>::scale(const T factor)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    for (const_it lambda = (*this).begin(); lambda != (*this).end(); ++lambda) {
        (*this).operator[]((*lambda).first) *= factor;
    }
}

template <typename T, typename Index>
void
Coefficients<Lexicographical,T,Index>::setToZero(void)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    for (const_it it = (*this).begin(); it!=(*this).end(); ++it) {
        (*this).operator[]((*it).first) = 0.;
    }
}

template <typename T, typename Index>
T
Coefficients<Lexicographical,T,Index>::norm(T tau) const
{
    //Coefficients<AbsoluteValue,T,Index> abs;
    //abs = *this;
    //return abs.norm(2.);
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    long double result=0.0L;
    if (Coefficients<Lexicographical,T,Index>::size() > 0) {
        for (const_it mu=Coefficients<Lexicographical,T,Index>::begin();
             mu!=Coefficients<Lexicographical,T,Index>::end(); ++mu) {
            result+=std::pow((long double)fabs((*mu).second), (long double)tau);
        }
    }
    return std::pow(result, (long double)(1.0L/tau));
}

template <typename T, typename Index>
std::ostream& operator<< (std::ostream &s, const Coefficients<Lexicographical,T,Index> &c)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    s << std::endl << "Coefficient_Rn<Lexicographical,T,Index>:" << std::endl;
    if (c.size() > 0) {
        for (const_it lambda = c.begin(); lambda != c.end(); ++lambda) {
            s << "  [" << (*lambda).first << "\t | " << (*lambda).second
                << "]" << std::endl;
        }
    }
    return s << std::endl;
}

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>
operator*(T alpha, const Coefficients<Lexicographical,T,Index> &_coeff) {
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    Coefficients<Lexicographical,T,Index> ret;
    for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
        ret.operator[]((*lambda).first) = alpha*(*lambda).second;
    }
    return ret;
}


template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>
P(const Coefficients<Lexicographical,T,Index> &v, const IndexSet<Index> &Lambda)
{
    typedef typename IndexSet<Index>::const_iterator const_set_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    Coefficients<Lexicographical,T,Index> ret;
    const_it v_end = v.end();
    for (const_set_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
        const_it it = v.find(*lambda);
        if (it != v_end) {
            ret[*lambda] = (*it).second;
        }
        else {
            ret[*lambda] = 0.;
        }
    }
    return ret;
}


template <typename T, typename Index>
void
P(const IndexSet<Index> &Lambda, Coefficients<Lexicographical,T,Index> &v)
{
    typedef typename IndexSet<Index>::const_iterator const_set_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_coeff_it;
    IndexSet<Index> indicesToBeDeleted(v.size());
    for (const_coeff_it it = v.begin(); it != v.end(); ++it) {
        if (Lambda.find((*it).first) == Lambda.end()) {
            indicesToBeDeleted.insert((*it).first);
        }
    }
    for (const_set_it lambda = indicesToBeDeleted.begin(); lambda != indicesToBeDeleted.end(); ++lambda) {
        v.erase(*lambda);
    }
    return;
}

template <typename T, typename Index>
IndexSet<Index>
supp(const Coefficients<Lexicographical,T,Index> &v)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    IndexSet<Index> ret;
    if (v.size() > 0) {
        for (const_it lambda = v.begin(); lambda != v.end(); ++lambda) {
            ret.insert((*lambda).first);
        }
    }

    return ret;
}

template <typename T, typename Index>
void
FillWithZeros(const IndexSet<Index> &Lambda, Coefficients<Lexicographical,T,Index> &u)
{
    typedef typename IndexSet<Index>::const_iterator const_it;
    for (const_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
        if (u.count((*lambda)) == 0) u[*lambda] = 0.0;
    }
}

template <typename T>
void
getLevelInfo(const Coefficients<Lexicographical,T,Index2D> &v, Index2D &maxIndex,
             Index2D &maxWaveletIndex, int *jmax, int &arrayLength)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    delete [] jmax;
    arrayLength = 2;
    jmax = new int[2];
    jmax[0] = -100; jmax[1] = -100;
    for (const_coeff_it it=v.begin(); it!=v.end(); ++it) {
        jmax[0] = std::max((int)(*it).first.index1.j,jmax[0]);
        jmax[1] = std::max((int)(*it).first.index2.j,jmax[1]);
        if ((*it).first.levelSum()>maxIndex.levelSum()) maxIndex = (*it).first;
        if ((*it).first.index1.xtype==XWavelet && (*it).first.index2.xtype==XWavelet) {
            if ((*it).first.levelSum()>maxWaveletIndex.levelSum()) maxWaveletIndex = (*it).first;
        }
        }
    return;
}

template <typename T>
void
getLevelInfo(const Coefficients<Lexicographical,T,Index3D> &v, Index3D &maxIndex,
             Index3D &maxWaveletIndex, int *jmax, int &arrayLength)
{
    typedef typename Coefficients<Lexicographical,T,Index3D>::const_iterator const_coeff_it;
    delete [] jmax;
    arrayLength = 3;
    jmax = new int[3];
    jmax[0] = -100; jmax[1] = -100; jmax[2] = -100;
    for (const_coeff_it it=v.begin(); it!=v.end(); ++it) {
        jmax[0] = std::max((int)(*it).first.index1.j,jmax[0]);
        jmax[1] = std::max((int)(*it).first.index2.j,jmax[1]);
        jmax[2] = std::max((int)(*it).first.index3.j,jmax[2]);
        if ((*it).first.levelSum()>maxIndex.levelSum()) maxIndex = (*it).first;
        if (    (*it).first.index1.xtype==XWavelet && (*it).first.index2.xtype==XWavelet
             && (*it).first.index3.xtype==XWavelet) {
            if ((*it).first.levelSum()>maxWaveletIndex.levelSum()) maxWaveletIndex = (*it).first;
        }
    }
}

template <typename T, typename Index>
void
Pe1(const Coefficients<Lexicographical,T,Index> &v,
          Coefficients<Lexicographical,T,Index1D> &Pe1_v)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator   const_coeff_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator const_coeff1d_it;
    for (const_coeff_it it=v.begin(); it!=v.end(); ++it) {
        const_coeff1d_it it_Pe1 = Pe1_v.find((*it).first.index1);
        if (it_Pe1==Pe1_v.end()) Pe1_v[(*it).first.index1] = 0.;
    }
}

template <typename T, typename Index>
void
Pe2(const Coefficients<Lexicographical,T,Index> &v,
          Coefficients<Lexicographical,T,Index1D> &Pe2_v)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator   const_coeff_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator const_coeff1d_it;
    for (const_coeff_it it=v.begin(); it!=v.end(); ++it) {
        const_coeff1d_it it_Pe2 = Pe2_v.find((*it).first.index2);
        if (it_Pe2==Pe2_v.end()) Pe2_v[(*it).first.index2] = 0.;
    }
}

template <typename T, typename Index>
void
Pe3(const Coefficients<Lexicographical,T,Index> &v,
          Coefficients<Lexicographical,T,Index1D> &Pe3_v)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator   const_coeff_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator const_coeff1d_it;
    for (const_coeff_it it=v.begin(); it!=v.end(); ++it) {
        const_coeff1d_it it_Pe3 = Pe3_v.find((*it).first.index3);
        if (it_Pe3==Pe3_v.end()) Pe3_v[(*it).first.index3] = 0.;
    }
}

/* **********************************************************************************************
 * Bucket sorted coefficient vector
 * ********************************************************************************************** */

template <typename T, typename Index>
Coefficients<Bucket,T,Index>::Coefficients()
: supremumnorm(0.), buckets(), bucket_ell2norms()
{

}

template <typename T, typename Index>
void
Coefficients<Bucket,T,Index>::bucketsort(const Coefficients<Lexicographical,T,Index> &_coeff, T eps)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    typedef typename Coefficients<Lexicographical,T,Index>::value_type val_type;

    for (int i=0; i<(int)buckets.size(); ++i) {
        buckets[i].clear();
    }
    buckets.clear();
    bucket_ell2norms.clear();

    supremumnorm = 0.;
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            supremumnorm = std::max(supremumnorm,(T)fabs((*lambda).second));
        }
    }
    //std::cerr << "Supremum norm = " << supremumnorm << std::endl;
    int NumOfBuckets = std::max(0,(int)(2*std::log(supremumnorm*std::sqrt(_coeff.size())/eps)/std::log(T(2)))+1);

    //std::cerr << "   Bucketsort: " << std::pow(2.,-0.5*NumOfBuckets)*supremumnorm*std::sqrt(_coeff.size())
    //          << " " << eps << std::endl;


    for (int i=0; i<NumOfBuckets; ++i) {
        BucketEntry tmp;
        buckets.push_back(tmp);
        bucket_ell2norms.push_back(0.);
    }

    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            int pos = std::max(0,int(std::ceil(-2*std::log(fabs((*lambda).second)/supremumnorm)/std::log(T(2))))-1);
            if (pos>=NumOfBuckets-1) continue;
            //std::cerr << "pos for " << (*lambda).first << ", " << (*lambda).second
            //          << ": " << -2*std::log(fabs((*lambda).second)/supremumnorm)/std::log(T(2)) << std::endl;
            const std::pair<const Index,T>* tmp = &(*lambda);
            //std::cout << "sizeof(lambda): " << sizeof(*lambda) << std::endl;
            //std::cout << "sizeof(tmo):    " << sizeof(tmp) << std::endl;

            buckets[pos].push_back(tmp);
            T val = (*lambda).second;
            bucket_ell2norms[pos] += (long double)(val*val);
        }
    }

    for (int i=0; i<(int)bucket_ell2norms.size(); ++i) {
        bucket_ell2norms[i] = std::sqrt(bucket_ell2norms[i]);
    }
}

template <typename T, typename Index>
int
Coefficients<Bucket,T,Index>::addBucketToIndexSet(IndexSet<Index> &Lambda, int bucketnumber)
{
    typedef typename  Coefficients<Bucket,T,Index>::BucketEntry::const_iterator const_it;
    for (const_it it=buckets[bucketnumber].begin(); it!=buckets[bucketnumber].end(); ++it) {
        //Index tmp((**it).first);
        Lambda.insert((**it).first);
    }
    return buckets[bucketnumber].size();
}

template <typename T, typename Index>
void
Coefficients<Bucket,T,Index>::addBucketToCoefficients(Coefficients<Lexicographical,T,Index> &coeff,
                                                      int bucketnumber)
{

    typedef typename  Coefficients<Bucket,T,Index>::BucketEntry::const_iterator const_it;
    for (const_it it=buckets[bucketnumber].begin(); it!=buckets[bucketnumber].end(); ++it) {
        Index tmp((**it).first);
        coeff[tmp] = (**it).second;
    }
}

template <typename T, typename Index>
std::ostream& operator<< (std::ostream &s, const Coefficients<Bucket,T,Index> &c)
{
    typedef typename  Coefficients<Bucket,T,Index>::BucketEntry::const_iterator const_it;
    s << std::endl << "Coefficients<Bucket,T,Index>:" << std::endl;
    if (c.buckets.size() > 0) {
        for (int i=0; i<(int)c.buckets.size(); ++i) {
            s << "Bucketnumber " << i+1 << " contains values in (" << c.supremumnorm*pow2ih<T>(-(i+1))
              << ", " << c.supremumnorm*pow2ih<T>(-i) << "]"
              << ", ell2norm=" << c.bucket_ell2norms[i] << std::endl;
            for (const_it it=c.buckets[i].begin(); it!=c.buckets[i].end(); ++it) {
                s << " " << (**it).first << " " << (**it).second << std::endl;
            }
        }
    }
    return s << std::endl;

}


/* **********************************************************************************************
 * Coefficient vector sorted by absolute values
 * ********************************************************************************************** */

template <typename T, typename Index>
Coefficients<AbsoluteValue,T,Index>::Coefficients()
{

}

template <typename T, typename Index>
Coefficients<AbsoluteValue,T,Index>&
Coefficients<AbsoluteValue,T,Index>::operator= (const Coefficients<Lexicographical,T,Index> &_coeff)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_it;
    typedef typename Coefficients<AbsoluteValue,T,Index>::value_type val_type;

    erase(Coefficients<AbsoluteValue,T,Index>::begin(), Coefficients<AbsoluteValue,T,Index>::end());
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            insert(val_type((*lambda).second, (*lambda).first));
        }
    }
    return *this;
}

template <typename T, typename Index>
Coefficients<AbsoluteValue,T,Index>&
Coefficients<AbsoluteValue,T,Index>::operator= (const Coefficients<AbsoluteValue,T,Index> &_coeff)
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_it;
    typedef typename Coefficients<AbsoluteValue,T,Index>::value_type val_type;

    erase(Coefficients<AbsoluteValue,T,Index>::begin(), Coefficients<AbsoluteValue,T,Index>::end());
    if (_coeff.size() > 0) {
        for (const_it lambda = _coeff.begin(); lambda != _coeff.end(); ++lambda) {
            insert(val_type((*lambda).first, (*lambda).second));
        }
    }
    return *this;
}

template <typename T, typename Index>
T
Coefficients<AbsoluteValue,T,Index>::norm(T tau) const
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_reverse_iterator const_it;
    T result=0.0;
    if (Coefficients<AbsoluteValue,T,Index>::size() > 0) {
        for (const_it mu=Coefficients<AbsoluteValue,T,Index>::rbegin();
             mu!=Coefficients<AbsoluteValue,T,Index>::rend(); ++mu) {
            result+=std::pow((T)fabs((*mu).first),(T) tau);
        }
    }

    return std::pow(result, 1.0/tau);
}

template <typename T, typename Index>
T
Coefficients<AbsoluteValue,T,Index>::l2bestnterm(int n) const
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_reverse_iterator const_rev_it;
    T result=0.0;
    int count=(*this).size();
    for (const_rev_it lambda=Coefficients<AbsoluteValue,T,Index>::rbegin(); lambda!=Coefficients<AbsoluteValue,T,Index>::rend(); ++lambda) {
        if (count>=n) {
            result += fabs((*lambda).first) * fabs((*lambda).first);
        }
        --count;
    }
    return sqrt(result);
}

template <typename T, typename Index>
T
Coefficients<AbsoluteValue,T,Index>::wtauNorm(T tau) const
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_it;
    T result=0.0, temp=0.0;
    int n=1;
    for (const_it lambda=Coefficients<AbsoluteValue,T,Index>::begin(); lambda!=Coefficients<AbsoluteValue,T,Index>::end(); ++lambda) {
        temp=std::pow(T(n),1.0/tau) * fabs((*lambda).first);
        if (temp>result) {
            result=temp;
        }
        ++n;
    }

    return result + norm();
}

template <typename T, typename Index>
DenseVector<Array<T> >
Coefficients<AbsoluteValue,T,Index>::norm_sections() const
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_it;

    DenseVector<Array<T> > result;

    if (Coefficients<AbsoluteValue,T,Index>::size() > 0) {
        //result.engine().resizeOrClear(_(0, int(log(T(Coefficients<AbsoluteValue,T,Index>::size()))/log(T(2))+1)));
        result.engine().resize(int(log(T(Coefficients<AbsoluteValue,T,Index>::size()))/log(T(2))+1));
        T help=0.0;
        int count=0, s=0;

        for (const_it lambda = Coefficients<AbsoluteValue,T,Index>::begin(); lambda != Coefficients<AbsoluteValue,T,Index>::end(); ++lambda) {
            help += ((*lambda).first)*((*lambda).first);
            count++;

            if (count==(1<<s)) {
                result(s+1) = sqrt(help);
                s++;
                help=0.0;
            }
        }
        if (help != 0.0) {
            result(s)=sqrt(help);
        }
    } else {
        result.engine().resize(0);
    }
    return result;
}

template <typename T, typename Index>
std::ostream& operator<< (std::ostream &s, const Coefficients<AbsoluteValue,T,Index> &c)
{
    typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator const_it;
    s << std::endl << "Coefficients<AbsoluteValue,T,Index>:" << std::endl;
    if (c.size() > 0) {
        for (const_it it=c.begin(); it!=c.end(); ++it) {
            s << "  [" << (*it).first << "\t | " << (*it).second << "]" << std::endl;
        }
    }
    return s << std::endl;
}


}   //namespace lawa

