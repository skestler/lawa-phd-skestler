#include <string.h>

namespace lawa {

template <typename T, typename RHSINTEGRAL, typename Preconditioner>
RHS1D<T,RHSINTEGRAL,Preconditioner>::RHS1D(const RHSINTEGRAL &_rhsintegral, Preconditioner &_P)
    :    rhsintegral(_rhsintegral), P(_P), rhs_data(), rhs_indexsets(), norm_estimate(0.),
         current_tol(2.), current_ell2norm(0.)
{

}

template <typename T, typename RHSINTEGRAL, typename Preconditioner>
bool
RHS1D<T,RHSINTEGRAL,Preconditioner>::readIndexSets(const char *filename)
{
    std::ifstream infile (filename);
    if (infile.is_open()) {
        std::cout << "File is open, ready to read..." << std::endl;

        std::string firstline, section, line;

        if (std::getline( infile, firstline, '\n' ) ) {
            norm_estimate = atof(firstline.c_str());
        }
        else {
            std::cerr << "Could not read first line" << std::endl;
            exit(1); return false;
        }

        while(std::getline( infile, section, '\n' )) {
            std::istringstream section_ss(section);
            std::string field;
            std::getline( section_ss, field, ',' );
            IndexSet<Index1D> Lambda;
            if (strcmp(field.c_str(),"#")==0) {
                std::getline( section_ss, field, ',' );
                std::cerr << "Reading index set for target accuray " << field << std::endl;

                std::string field1, field2, field3;
                while(std::getline( infile, line, '\n' )) {
                    std::istringstream line_ss(line);
                    std::getline( line_ss, field1, ',' );
                    std::getline( line_ss, field2, ',' );
                    std::getline( line_ss, field3, ',' );
                    int j,k;

                    if (strcmp(field1.c_str(),"wavelet")==0) {
                        j = atoi(field2.c_str());
                        k = atoi(field3.c_str());
                        Lambda.insert(Index1D(j,k,XWavelet));
                        //std::cout << "wavelet, " << field2 << ", " << field3 << std::endl;
                    }
                    else if (strcmp(field1.c_str(),"scaling")==0) {
                        j = atoi(field2.c_str());
                        k = atoi(field3.c_str());
                        Lambda.insert(Index1D(j,k,XBSpline));
                        //std::cout << "scaling, " << field2 << ", " << field3 << std::endl;
                    }
                    else if (strcmp(field1.c_str(),"")==0) {
                        rhs_indexsets.push_back(Lambda);
                        break;
                    }
                    else {
                        std::cerr << "Got " << field1 << ", could not read file." << std::endl;
                        exit(1); return false;
                    }
                }
            }
            else {
                std::cerr << "Got " << field << ", could not read file." << std::endl;
                exit(1); return false;
            }
        }
        return true;
    }
    else {
        return false;
    }

}

template <typename T, typename RHSINTEGRAL, typename Preconditioner>
T
RHS1D<T,RHSINTEGRAL,Preconditioner>::operator()(const Index1D &lambda)
{
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator const_coeff_it;
    typedef typename Coefficients<AbsoluteValue,T,Index1D>::value_type val_type;
    const_coeff_it it_end       = rhs_data.end();
    const_coeff_it it_index     = rhs_data.find(lambda);

    if (it_index != it_end) {
        return (*it_index).second;
    }
    else {
        T ret = P(lambda) * rhsintegral(lambda);
        rhs_data[lambda] = ret;
        return ret;
    }
}

template <typename T, typename RHSINTEGRAL, typename Preconditioner>
Coefficients<Lexicographical,T,Index1D>
RHS1D<T,RHSINTEGRAL,Preconditioner>::operator()(const IndexSet<Index1D> &Lambda)
{
    typedef typename IndexSet<Index1D>::const_iterator const_set_it;
    Coefficients<Lexicographical,T,Index1D> ret;
    for (const_set_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
        T tmp = this->operator()(*lambda);
        ret[*lambda] = tmp;
    }
    return ret;
}

template <typename T, typename RHSINTEGRAL, typename Preconditioner>
Coefficients<Lexicographical,T,Index1D>
RHS1D<T,RHSINTEGRAL,Preconditioner>::operator()(T tol)
{
    if (tol >= current_tol) {
        Coefficients<Lexicographical,T,Index1D> ret;
        Coefficients<Bucket,T,Index1D>          rhs_bucket;

        T thresh = (tol-current_tol)/std::sqrt(rhs_data.size());
        //std::cerr << "current_tol = " << current_tol << ", thresh = " << thresh << std::endl;
        //std::cerr << "current_ell2norm = " << current_ell2norm << std::endl;

        rhs_bucket.bucketsort(rhs_data,thresh);

        //long double squared_ell2norm=current_ell2norm*current_ell2norm;
        long double squared_ell2norm=0.0L;
        for (int i=0; i<(int)rhs_bucket.bucket_ell2norms.size(); ++i) {
            squared_ell2norm += (long double)std::pow(rhs_bucket.bucket_ell2norms[i],2.L);
            rhs_bucket.addBucketToCoefficients(ret,i);
            //std::cerr << "(" << i << ", " << rhs_bucket.bucket_ell2norms.size() << ")" << std::endl;
            //std::cerr << "  -> ell2-norm: " << squared_ell2norm << std::endl;
            if (fabs(current_ell2norm*current_ell2norm-squared_ell2norm) <= std::pow(tol-current_tol,2.)) break;
        }
        return ret;

    }
    else {
        std::cerr << "RHS1D<...> : Decreasing current tolerance from " << current_tol << " to "
                  << tol << std::endl;
        long double squared_ell2norm=current_ell2norm*current_ell2norm;
        while(current_tol>tol) {
            if (rhs_indexsets.size()==0 && current_tol>tol) {
                std::cerr << "Attention: Insufficient rhs information." << std::endl;
                exit(1);
            }
            if (rhs_indexsets.size()>0) {
                for (const_set1d_it it = (*rhs_indexsets.begin()).begin();
                                    it!= (*rhs_indexsets.begin()).end(); ++it) {
                    T val = rhsintegral(*it)*P(*it);
                    rhs_data[(*it)] = val;
                    squared_ell2norm += (long double)val*val;
                }
                rhs_indexsets.pop_front();
            }
            else {
                continue;
            }
            current_tol *= 0.5;
        }
        current_ell2norm = std::sqrt(squared_ell2norm);

        return this->operator()(tol);
    }

}

}   // namespace lawa
