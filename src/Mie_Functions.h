

int ConfigReader(const char *filename, double D[]);
vector<vector<double>> RectPtsGen(double Lx, double Ly, double Lz, double Cx, double Cy, double Cz, double res);
vector<vector<double>> SpherePtsGen(double R, double res);
vector<vector<double>> CirclePtsGen(int plane, double R, double res);
complex<double> J1(double n, complex<double> z);
double Y1(double fnu, double x, double k);
complex<double> J1d(double n, complex<double> z);
double Y1d(double fnu, double x, double k);
complex<double> H1(double fnu, complex<double> z);
complex<double> H2(double n, complex<double> z);
complex<double> H2d(double n, complex<double> z);
double P0(double n, double theta);
double P1(double n, double theta);
double P1d(double n, double theta);
double CscaPEC(double N, double r, double lam0, complex<double> erb, complex<double> mrb);
complex<double> * Mie_fieldsPEC(double IN, double x, double y, double z, double N, double r, double lam0, complex<double> erb, complex<double> mrb, complex<double> Fields[]);
complex<double> * TM_fieldsPEC(double N, double R, double r, double k0, double theta, double phi, complex<double> TM0[]);
complex<double> * TE_fieldsPEC(double N, double R, double r, double k0, double theta, double phi, complex<double> TE0[]);
double * force(vector<vector<complex<double>>> TM, vector<vector<complex<double>>> TE, vector<vector<double>> angle, double dA, double R, double F[]);
void CsToFileR(const char *Name, vector<double> rv, vector<complex<double>> CS, double lam, double n, complex<double> erb, complex<double> mrb);
void CsToFileW(const char *Name, vector<double> wv, vector<complex<double>> CS, double r, double n, complex<double> erb, complex<double> mrb);
void FieldsToFileXYZ(const char *Name, vector<vector<double>> A, vector<vector<complex<double>>> vT, double r, double lam, double n, complex<double> erb, complex<double> mrb);
void FieldsToFile(const char *Name, vector<vector<double>> A, vector<vector<complex<double>>> vT);
void ForcesToFile(const char *Name, vector<vector<double>> F, vector<double> r, double D[]);
vector<double> range(double start, double stop, double step);
vector<double> linspace(double start, double stop, int N);


