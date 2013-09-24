namespace lawa {

template<typename T, typename Index>
void
writeCoefficientsToFile(Coefficients<Lexicographical,T,Index> &u, int i, const char* filename)
{
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator    const_coeff_it;

    std::stringstream filenamestr;
    filenamestr << filename << "__" << i << ".dat";
    std::ofstream file(filenamestr.str().c_str());
    file.precision(20);
    std::cerr << "Started writing into file..." << std::endl;
    for (const_coeff_it it=u.begin(); it!=u.end(); ++it) {
        file << (*it).first << " " << (*it).second << std::endl;
    }
    file.close();
    std::cerr << "... finished." << std::endl;
}


template<typename T>
void
readCoefficientsFromFile(Coefficients<Lexicographical,T,Index2D> &u, const char* filename)
{
    std::ifstream infile (filename);
    if (infile.is_open()) {
        std::cout << "File is open, ready to read..." << std::endl;
        std::string line;
        //Coefficients<Lexicographical,T,Index1D> tmp;
        int count = 0;
        while(std::getline( infile, line, '\n' )) {
            //cout << line << endl;
            std::string field1, field2, field3, field4, field5, field6, field7;
            std::istringstream line_ss(line);
            std::getline( line_ss, field1, ',' );
            std::getline( line_ss, field2, ',' );
            std::getline( line_ss, field3, ',' );
            std::getline( line_ss, field4, ',' );
            std::getline( line_ss, field5, ',' );
            std::getline( line_ss, field6, ' ' );
            std::getline( line_ss, field7, ' ' );
            //std::cerr << field1 << " " << field2 << " " << field3 << " " << field4
            //          << " " << field5 << " " << field6 << " " << field7 << std::endl;
            Index1D index1, index2;
            double val;
            int j1 = atoi(field2.c_str());
            int k1 = atoi(field3.c_str());
            index1.j = j1; index1.k = k1;
            int j2 = atoi(field5.c_str());
            int k2 = atoi(field6.c_str());
            index2.j = j2; index2.k = k2;
            val = atof(field7.c_str());
            //std::cerr << field1 << " " << j1 << " " << k1 << " " << field4 << " " << j2 << " " << k2 << " " << val << " " << field7 << endl << endl;

            if (strcmp(field1.c_str(),"wavelet")==0)    index1.xtype = XWavelet;
            else                                        index1.xtype = XBSpline;
            if (strcmp(field4.c_str(),"wavelet")==0)    index2.xtype = XWavelet;
            else                                        index2.xtype = XBSpline;

            Index2D index(index1,index2);
            u[index] = val;
            ++count;
        }
        std::cerr << "Read " << count << " coefficients, #supp u = " << u.size() << std::endl;
        //cout << "Size of 1d-tree: " << tmp.size() << endl;
        //plotCoeff(tmp, basis, "coeff", false, true);
    }
    else {
        std::cout << "File " << filename << " not found." << std::endl;
        exit(1);
        return;
    }
    //cout << "u = " << u << endl;
}

template<typename T>
void
readCoefficientsFromFile(Coefficients<Lexicographical,T,Index3D> &u, const char* filename)
{
    std::ifstream infile (filename);
    if (infile.is_open()) {
        std::cout << "File " << filename << " is open, ready to read..." << std::endl;
        std::string line;
        //Coefficients<Lexicographical,T,Index1D> tmp;
        int count = 0;
        while(std::getline( infile, line, '\n' )) {
            //cout << line << endl;
            std::string field1,field2,field3,field4,field5,field6,field7,field8,field9,field10;
            std::istringstream line_ss(line);
            std::getline( line_ss, field1, ',' );
            std::getline( line_ss, field2, ',' );
            std::getline( line_ss, field3, ',' );
            std::getline( line_ss, field4, ',' );
            std::getline( line_ss, field5, ',' );
            std::getline( line_ss, field6, ',' );
            std::getline( line_ss, field7, ',' );
            std::getline( line_ss, field8, ',' );
            std::getline( line_ss, field9, ' ' );
            std::getline( line_ss, field10, ' ' );
            //std::cerr << field1 << " " << field2 << " " << field3 << " "
            //          << field4 << " " << field5 << " " << field6 << " "
            //          << field7 << " " << field8 << " " << field9 << " " << field10 << std::endl;
            Index1D index1, index2, index3;
            double val;
            int j1 = atoi(field2.c_str());
            int k1 = atoi(field3.c_str());
            index1.j = j1; index1.k = k1;
            int j2 = atoi(field5.c_str());
            int k2 = atoi(field6.c_str());
            index2.j = j2; index2.k = k2;
            int j3 = atoi(field8.c_str());
            int k3 = atoi(field9.c_str());
            index3.j = j3; index3.k = k3;
            val = atof(field10.c_str());
            //std::cerr << field1 << " " << j1 << " " << k1 << " " << field4 << " " << j2 << " " << k2 << " " << val << " " << field7 << endl << endl;

            if (strcmp(field1.c_str(),"wavelet")==0)    index1.xtype = XWavelet;
            else                                        index1.xtype = XBSpline;
            if (strcmp(field4.c_str(),"wavelet")==0)    index2.xtype = XWavelet;
            else                                        index2.xtype = XBSpline;
            if (strcmp(field7.c_str(),"wavelet")==0)    index3.xtype = XWavelet;
            else                                        index3.xtype = XBSpline;
            ++count;
            Index3D index(index1,index2,index3);
            u[index] = val;
        }
        std::cerr << "   Read " << count << " coefficients, #supp u = " << u.size() << std::endl;
        //cout << "Size of 1d-tree: " << tmp.size() << endl;
        //plotCoeff(tmp, basis, "coeff", false, true);
    }
    else {
        std::cout << "File " << filename << " not found." << std::endl;
        exit(1);
        return;
    }
    infile.close();
    //cout << "u = " << u << endl;
}

}   // namespace lawa
