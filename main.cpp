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
	int new_n;
	float t_simple;
	float t_fft;
	vector <int> A;
	vector <int> B;
	vector <double> C;
public:
	fft(int &);
	void setupDegree(int &);
	void simpleMultiplication(int &);
	vector< complex<double> > DFT (vector <int> &);
	void interpolation (vector <complex <double> > &);
	void fastFourierTransformation(int &);

};

fft::fft(int &n)
{
	int k;				//temp. condition
	cout<<"Enter coeff. of A(x) (starting from lower order including 0): \n";
	for(int i = 0; i<n; i++)
	{
		cin>>k;
		A.push_back(k);
	}

	cout<<"Enter coeff. of B(x): \n";
	for(int i = 0; i<n; i++)
	{
		cin>>k;
		B.push_back(k);
	}

	setupDegree(n);
	return;
}

void fft::setupDegree(int &n)
{
	int i=0;
	int k = ((2*n)-1);
	while(1)
	{
		if(pow(2,i)>=k)
			break;
		i++;
	}
	new_n = pow(2, i);
	return;
}

void fft::simpleMultiplication(int &n)
{
	int temp=0;
	clock_t t;
	t = clock();
	for(int i = 0; i<(2*n - 1); i++)			//(2n-1) is the degree bound here!
	{
		for(int j = 0; j<=i; j++)
		{
			if(j<n && (i-j)<n)
				temp += A[j]*B[i-j];
		}
		C.push_back(temp);
		temp = 0;
	}
	t = clock() - t;
	t_simple = ((float)t)/CLOCKS_PER_SEC;
		#ifdef DEBUG
		cout<<"MULTI_SIMPLE"<<endl;
		for(int i = 0; i<(2*n - 1); i++)
			cout<<C[i]<<" ";
        cout<<endl;
		#endif
	C.resize(0);
	return;
}

vector < complex<double> > fft::DFT (vector <int> &V)
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


		vector < int > a1 ;
		vector < int > a2 ;

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

void fft::interpolation (vector <complex <double> > &CP)
{
	complex <double> omega_n (cos(2*PI/new_n),sin(2*PI/new_n));
	vector < complex <double> > temp;
	for (int i = 0; i < new_n; ++i)
	{
	    complex <double> k (0,0);
		for (int j = 0; j < new_n; ++j)
		{
			k += (CP[j]*(pow(omega_n,-(i*j))));
		}
		C.push_back(k.real()/new_n);
	}
}

void fft::fastFourierTransformation(int &n)
{
	clock_t t;
	t = clock();
	vector <complex <double> >Y1;
	vector <complex <double> >Y2;
	while(A.size() < new_n)
		A.push_back(0);
	while(B.size() < new_n)
		B.push_back(0);
    cout<<endl;
	Y1 = DFT(A);
	Y2 = DFT(B);

										//Multiplication [step 3]
	vector <complex <double> >C_Pointwise;
	for (int i = 0; i < new_n; ++i)
		C_Pointwise.push_back(Y1[i]*Y2[i]);

	interpolation(C_Pointwise);

	t = clock() - t;
	t_fft = ((float)t)/CLOCKS_PER_SEC;

	#ifdef
	cout<<"MULTI_FFT"<<endl;
	for (int i=0; i<(2*n - 1); i++)
        cout<<C[i]<<" ";
    cout<<endl;
    #endif
    cout<<t_simple<<" "<<t_fft<<endl;



	#ifdef DEBUG
	cout<<"check\n";
	for (int i = 0; i < Y1.size(); ++i)
		cout<<Y1[i]<<" "<<Y2[i]<<endl;
	#endif

}



int main()
{
	int n;
	cout<<"degree bound?\n";
		cin>>n;
	fft f(n);
	f.simpleMultiplication(n);
	f.fastFourierTransformation(n);
	return 0;
}
