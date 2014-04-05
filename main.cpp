#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <ctime>
#define PI 3.1415926535897
using namespace std;

class fft
{
private:
	int na,nb;
	int new_n;
	float t_fft;
	vector <double> A;
	vector <double> B;
	vector <double> C;
	vector <complex <double> >C_Pointwise;
public:
	fft(int &, int &);
	void setupDegree();
	void finalise();
	vector< complex<double> > DFT (vector <double> &);
	vector< complex<double> > interpolation (vector< complex<double> > &);
	void fastFourierTransformation();

};

fft::fft(int &n1, int &n2)
{
	na=n1;
	nb=n2;
	int k;				//temp. condition
	cout<<"Enter coeff. of A(x) (starting from lower order including 0): \n";
	for(int i = 0; i<na; i++)
	{
		cin>>k;
		A.push_back(k);
	}

	cout<<"Enter coeff. of B(x): \n";
	for(int i = 0; i<nb; i++)
	{
		cin>>k;
		B.push_back(k);
	}
	return;
}

void fft::setupDegree()
{
	int i=0;
	int k = (na+nb-1);
	while(1)
	{
		if(pow(2,i)>=k)
			break;
		i++;
	}
	new_n = pow(2, i);
	return;
}

vector < complex<double> > fft::DFT (vector <double> &V)
{
	int n = V.size();

	if(n==1)
	{
		vector <complex <double> >R;
		complex <double> temp (V[0],0);
		R.push_back(temp);
		return R;
	}

	else
    {
        complex <double> omega_n (cos(2*PI/n),sin(2*PI/n));			//w_n
		complex <double> omega (1,0);								//w=1

		vector < complex<double> > v1 ;
		vector < complex<double> > v2 ;
		vector < complex<double> > y (n) ;


		vector < double > a1 ;
		vector < double > a2 ;

		for (int i = 0; i < n; i+=2)
		{
			a1.push_back(V[i]);
			a2.push_back(V[i+1]);
		}

		v1 = DFT(a1);
		v2 = DFT(a2);

		for (int i = 0; i < (n/2); ++i)
		{
			y[i] = v1[i] + omega*v2[i];
			y[i+(n/2)] = v1[i] - omega*v2[i];
			omega *= omega_n;
		}
		return y;
    }
}


vector< complex<double> > fft::interpolation (vector< complex<double> > &V)
{
	int n = V.size();

	if(n==1)
	{
		vector <complex <double> >R;
		complex <double> temp (real(V[0]),imag(V[0]));
		R.push_back(temp);
		return R;
	}

	else
    {
      complex <double> omega_n (cos(2*PI/n),-sin(2*PI/n));				//w_n
		complex <double> omega (1,0);												//w=1

		vector < complex<double> > v1 ;
		vector < complex<double> > v2 ;
		vector < complex<double> > y (n) ;


		vector < complex<double> > a1 ;
		vector < complex<double> > a2 ;

		for (int i = 0; i < n; i+=2)
		{
			a1.push_back(V[i]);
			a2.push_back(V[i+1]);
		}

		v1 = interpolation(a1);
		v2 = interpolation(a2);

		for (int i = 0; i < (n/2); ++i)
		{
			y[i] = (v1[i] + omega*v2[i]);
			y[i+(n/2)] = (v1[i] - omega*v2[i]);
			omega *= omega_n;
		}
		return y;
    }
}

void fft::finalise()
{
	for (int i=0; i<(na + nb - 1); i++)
		C.push_back ( real(C_Pointwise[i] ) / new_n);
	return;
}

void fft::fastFourierTransformation()
{
	clock_t t;
	t = clock();
	setupDegree();
	vector <complex <double> >Y1;
	vector <complex <double> >Y2;

	//Adding zeros for rest of the term till it becomes closest 2^n
	while(A.size() < new_n)
		A.push_back(0);
	while(B.size() < new_n)
		B.push_back(0);

	//Evaluation
	Y1 = DFT(A);
	Y2 = DFT(B);
	
	
	//Point wise multiplication
	for (int i = 0; i < new_n; ++i)
		C_Pointwise.push_back(Y1[i]*Y2[i]);

	//Interpolation
	C_Pointwise = interpolation(C_Pointwise);
	finalise();

	t = clock() - t;
	t_fft = ((float)t)/CLOCKS_PER_SEC;

	cout<<"MULTI_FFT"<<endl;
	for (int i=0; i<(na + nb - 1); i++)
        cout<<C[i]<<" ";
    cout<<endl;
    cout<<"Time: "<<t_fft<<endl;

	#ifdef DEBUG
	cout<<"check\n";
	for (int i = 0; i < Y1.size(); ++i)
		cout<<Y1[i]<<" "<<Y2[i]<<endl;
	#endif
}

int main()
{
	int na,nb;
	cout<<"degree bound of A and B respectively?\n";
		cin>>na>>nb;
	fft f(na, nb);
	f.fastFourierTransformation();
	return 0;
}