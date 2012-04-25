#ifndef __INCLUDE_MEASURED_DATA
#define __INCLUDE_MEASURED_DATA

#include <iomanip>

std::vector<ComplexType> data1GM00(n_zone*wn_max*n_part,0);
extern const int n_part;
extern const int n_zone;
extern int wn_max;

#ifdef use_mpi
    #define output_stream CTQMC.getStream()
#else
    #define output_stream std::cout
#endif

struct MeasuredData
{
    std::vector<n_type> dataGM00_real; //n_zone*wn_max*n_part);//,0);
    std::vector<n_type> dataGM00_imag; //n_zone*wn_max*n_part);//,0);
    std::vector<n_type> dataGtotal_real;// (n_zone*wn_max*n_part*n_part,0);
    std::vector<n_type> dataGtotal_imag;// (n_zone*wn_max*n_part*n_part,0);
    std::vector<n_type> dataGamma4_real;// (n_zone*wn_max*n_part*n_part,0);
    std::vector<n_type> dataGamma4_imag;// (n_zone*wn_max*n_part*n_part,0);
    n_type weight_sum;
    
    MeasuredData();

    void setFromComplexData(std::vector<ComplexType> &inGM00, std::vector<ComplexType> &inGtotal);
    void setGamma(std::vector<ComplexType> &inGamma);
    
    std::string& serializeGM00();
    std::string& serializeGtotal();
};

MeasuredData::MeasuredData()
{
    dataGM00_real.assign(n_zone*wn_max*n_part,0.0);
    dataGM00_imag.assign(n_zone*wn_max*n_part,0.0);
    dataGtotal_real.assign(n_zone*wn_max*n_part*n_part,0.0);
    dataGtotal_imag.assign(n_zone*wn_max*n_part*n_part,0.0);
    weight_sum=0;
    output_stream << "Initialized data containers" << std::endl;
};

void MeasuredData::setFromComplexData(std::vector<ComplexType> &inGM00, std::vector<ComplexType> &inGtotal)
{
for (int w=0; w<wn_max; ++w){
    for (int z=0; z<n_zone; ++z){
        for (int i=0; i<n_part; ++i){
            dataGM00_real[z*wn_max*n_part + w*n_part+i] = real(inGM00[z*wn_max*n_part + w*n_part+i]);
            dataGM00_imag[z*wn_max*n_part + w*n_part+i] = imag(inGM00[z*wn_max*n_part + w*n_part+i]);
            for (int j=0; j<n_part; ++j){
                dataGtotal_real[z*wn_max*n_part*n_part + w*n_part*n_part+i*n_part+j] = real(inGtotal[z*wn_max*n_part*n_part + w*n_part*n_part+i*n_part+j]);
                dataGtotal_imag[z*wn_max*n_part*n_part + w*n_part*n_part+i*n_part+j] = imag(inGtotal[z*wn_max*n_part*n_part + w*n_part*n_part+i*n_part+j]);
                }; // end of j loop
            }; // end i loop
        }; // end of z loop
    }; // end of w loop
}

std::string& MeasuredData::serializeGM00()
{
static std::stringstream output;
for (int w=0; w<wn_max; ++w){
    output << w << "   " << std::flush;
    for (int z=0; z<n_zone; ++z){
        for (int i=0; i<n_part; ++i){
            output << std::fixed << std::setprecision(9) << dataGM00_real[z*wn_max*n_part + w*n_part+i]/weight_sum << " " 
                                                         << dataGM00_imag[z*wn_max*n_part + w*n_part+i]/weight_sum << "  " << std::flush;
            }; // end i loop
        output << "  " << std::flush;
        }; // end of z loop
    output << std::endl;
    }; // end of w loop
static std::string out_str = output.str();
return out_str;
};

std::string& MeasuredData::serializeGtotal()
{
static std::stringstream output;
static std::string out_string;
for (int w=0; w<wn_max; ++w){
    output << w << "   " << std::flush;
    for (int z=0; z<n_zone; ++z){
        for (int i=0; i<n_part; ++i){
            for (int j=0; j<n_part; ++j){
                output << std::fixed << std::setprecision(9)
                    << dataGtotal_real[z*wn_max*n_part*n_part + w*n_part*n_part+i*n_part+j]/weight_sum << " " 
                    << dataGtotal_imag[z*wn_max*n_part*n_part + w*n_part*n_part+i*n_part+j]/weight_sum <<  "  " << std::flush;
                }; // end of j loop
            }; // end i loop
            output << "  " << std::flush;
        }; // end of z loop
    output << std::endl;
    }; // end of w loop
out_string=output.str();
return out_string;
};

#endif // end::#ifndef __INCLUDE_MEASURED_DATA
