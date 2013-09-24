#include <string.h>

namespace lawa {

template <typename T, typename RHSINTEGRAL, typename Preconditioner>
RHS2D<T,RHSINTEGRAL,Preconditioner>::RHS2D(const RHSINTEGRAL &_rhsintegral, Preconditioner &_P)
    :    rhsintegral(_rhsintegral), P(_P), rhs_data(), rhs_indexsets(), norm_estimate(0.),
         current_tol(2.), current_ell2norm(0.), new_values(false)
{

}

template <typename T, typename RHSINTEGRAL, typename Preconditioner>
bool
RHS2D<T,RHSINTEGRAL,Preconditioner>::readIndexSets(const char *filename)
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
            IndexSet<Index2D> Lambda;
            if (strcmp(field.c_str(),"#")==0) {
                std::getline( section_ss, field, ',' );
                std::cerr << "Reading index set for target accuray " << field << std::endl;

                std::string field1, field2, field3, field4, field5, field6;
                while(std::getline( infile, line, '\n' )) {
                    std::istringstream line_ss(line);
                    std::getline( line_ss, field1, ',' );
                    std::getline( line_ss, field2, ',' );
                    std::getline( line_ss, field3, ',' );
                    std::getline( line_ss, field4, ',' );
                    std::getline( line_ss, field5, ',' );
                    std::getline( line_ss, field6, ',' );
                    int j1,j2;
                    long k1,k2;

                    j1 = atoi(field2.c_str());
                    k1 = atol(field3.c_str());
                    j2 = atoi(field5.c_str());
                    k2 = atol(field6.c_str());

                    if (strcmp(field1.c_str(),"wavelet")==0 && strcmp(field4.c_str(),"wavelet")==0) {
                        Index1D index_x(j1,k1,XWavelet);
                        Index1D index_y(j2,k2,XWavelet);
                        Lambda.insert(Index2D(index_x,index_y));
                    }
                    else if (strcmp(field1.c_str(),"wavelet")==0 && strcmp(field4.c_str(),"scaling")==0) {
                        Index1D index_x(j1,k1,XWavelet);
                        Index1D index_y(j2,k2,XBSpline);
                        Lambda.insert(Index2D(index_x,index_y));
                    }
                    else if (strcmp(field1.c_str(),"scaling")==0 && strcmp(field4.c_str(),"wavelet")==0) {
                        Index1D index_x(j1,k1,XBSpline);
                        Index1D index_y(j2,k2,XWavelet);
                        Lambda.insert(Index2D(index_x,index_y));
                    }
                    else if (strcmp(field1.c_str(),"scaling")==0 && strcmp(field4.c_str(),"scaling")==0) {
                        Index1D index_x(j1,k1,XBSpline);
                        Index1D index_y(j2,k2,XBSpline);
                        Lambda.insert(Index2D(index_x,index_y));
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
        std::cerr << "Could not read file " << filename << std::endl;
        exit(1);
        return false;
    }

}

template <typename T, typename RHSINTEGRAL, typename Preconditioner>
T
RHS2D<T,RHSINTEGRAL,Preconditioner>::operator()(const Index2D &lambda)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    typedef typename Coefficients<AbsoluteValue,T,Index2D>::value_type val_type;
    const_coeff_it it_end       = rhs_data.end();
    const_coeff_it it_index     = rhs_data.find(lambda);

    if (it_index != it_end) {
        return (*it_index).second;
    }
    else {
        T ret = P(lambda) * rhsintegral(lambda);
        rhs_data[lambda] = ret;
        new_values = true;
        return ret;
    }
}

template <typename T, typename RHSINTEGRAL, typename Preconditioner>
Coefficients<Lexicographical,T,Index2D>
RHS2D<T,RHSINTEGRAL,Preconditioner>::operator()(const IndexSet<Index2D> &Lambda)
{
    typedef typename IndexSet<Index2D>::const_iterator const_set_it;
    Coefficients<Lexicographical,T,Index2D> ret;
    for (const_set_it lambda = Lambda.begin(); lambda != Lambda.end(); ++lambda) {
        T tmp = this->operator()(*lambda);
        ret[*lambda] = tmp;
    }
    return ret;
}

template <typename T, typename RHSINTEGRAL, typename Preconditioner>
Coefficients<Lexicographical,T,Index2D>
RHS2D<T,RHSINTEGRAL,Preconditioner>::operator()(T tol)
{
    if (tol > current_tol  || new_values) {
        new_values = false;
        Coefficients<Lexicographical,T,Index2D> ret;
        Coefficients<Bucket,T,Index2D>          rhs_bucket;

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
            if (fabs(current_ell2norm*current_ell2norm-squared_ell2norm) <= std::pow(tol-current_tol,(T)2.)) break;
        }
        return ret;

    }
    else {
        std::cerr << "RHS2D<...> : Decreasing current tolerance from " << current_tol << " to "
                  << tol << std::endl;
        long double squared_ell2norm=current_ell2norm*current_ell2norm;
        while(current_tol>=tol) {
            if (rhs_indexsets.size()>0) {
                for (const_set_it it = (*rhs_indexsets.begin()).begin();
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
        if (rhs_indexsets.size()==0 && current_tol>=tol) {
            std::cerr << "Attention: Insufficient rhs information." << std::endl;
            exit(1);
        }
        return this->operator()(tol);
    }
}

template <typename T, typename RHSINTEGRAL, typename Preconditioner>
IndexSet<Index2D>
RHS2D<T,RHSINTEGRAL,Preconditioner>::getFullIndexSet()
{
    typedef std::list<IndexSet<Index2D> >::const_iterator const_list_it;
    IndexSet<Index2D> Lambda;
    for (const_list_it it_list=rhs_indexsets.begin(); it_list!=rhs_indexsets.end(); ++it_list) {
        std::cout << " running through indexset list" << std::endl;
        for (const_set_it it = (*it_list).begin(); it!= (*it_list).end(); ++it) {
            Lambda.insert(*it);
        }
    }
    return Lambda;
}

}   // namespace lawa
