#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <math.h>
#include <complex>
#include <iomanip>

using namespace std;

void meanCentre(complex <float> **datapoint, const int count, const int leads);
void t(int i);
void writeWTplt(int leads, string patient);
void writeFTplt(int leads, string patient);
void waveletRange(float p, int k, int timebins, int *mMin, int *mMax);
void chooseWavelet(string *waveletChoice, const float sds);
void chooseTransform(string *transformChoice);
/*WAVELET PROTOTYPES*/
float gauss(complex <float> pi, float x, float mu, float sigma);
float ricker(complex <float> pi, float x, float mu, float sigma);
float lorentz(complex <float> pi, float x, float mu, float sigma);
float haar(const int i, const float p);
/*END WAVELET PROTOTYPES*/
void determineBins(int count, int n, int * timebins, int * n_lastbin);

class spectrum
{
	protected:
		float FREQUENCY;
		complex <float> FOURIER;
	public:
		spectrum()
		{
			FREQUENCY = .0;
			FOURIER.real() = .0;
			FOURIER.imag() = .0;
		}
		spectrum(float frequency, complex <float> fourier)
		{
			FREQUENCY = frequency;
			FOURIER = fourier;
		}
		complex <float> FT()
		{return FOURIER;}

		float w()
		{return FREQUENCY;}

		float a2()
		{return (FOURIER*conj(FOURIER)).real();}

		spectrum operator+(const complex <float> &x) const //spectrum(freq,FT) + complex <float> => spectrum(freq,FT+complex <float>)
		{
			spectrum temp(FREQUENCY,FOURIER+x);
			return temp;
		}

		~spectrum(){}
};

int main()
{
	const int leads(12);
	int patients(0);
	int * patientArray;
	float **spreadArray;

	cout << "12-Lead ECG Fourier/Wavelet Transform" << endl;
	cout << "Enter number of patients:" << endl;
	//cin >> patients;
	patients = 1;
	patientArray = new int [patients];
	for(int x(0); x < patients; x = x + 1)
	{
		cout << "Enter patient (patientArray[" << x << "]) index:" << endl;
		cin >> patientArray[x];
	}

	spreadArray = new float * [patients];
	for(int patient = 0; patient < patients; patient = patient + 1)
	{
		spreadArray[patient] = new float [leads]();
	}
for(int z(0); z < patients; z = z + 1)
{
	/*INT DECLARATIONS*/
	int patient(0);
	int count(0);
	int sampleCount(0);
	int n(0);
	int n_lastbin(0);
	int timebins(0);
	int freqbins(0);
	int zeroPads(0);
	/*END INT DECLARATIONS*/

	/*FLOAT DECLARATIONS*/
	float p(.0);
	float dt(.0);	
	float datum(.0);
	float arg(.0);	
	float w_i(.0);
	float w_f(.0);
	float dw(.0);
	float f_i(.0);
	float f_f(.0);
	float df(.0);
	float sigma(.0);
	float mu(.0);
	float sds(4.5);
	float nyquist(0);
	float sumFT(0);
	float sumFT2(0);
	float spread(0);
	/*END FLOAT DECLARATIONS*/

	/*STRING DECLARATIONS*/
	string transformChoice("");
	string waveletChoice("");
	string sourceDataFilename("");
	string input("");
	string denoiseChoice("");
	string plotChoice("");
	/*END STRING DECLARATIONS*/

	/*BOOL DECLARATIONS*/
	bool binary(1);
	/*END BOOL DECLARATIONS*/

	/*COMPLEX <FLOAT> DECLARATIONS*/
	const complex <float> pi(3.14159,.0);
	complex <float> cbins(.0,.0);
	complex <float> sum(.0,.0);
	complex <float> mean(.0,.0);
	complex <float> expo(.0,.0);
	complex <float> ccount(.0,.0);
	complex <float> clocalFunction(.0,.0);
	complex <float> cn(.0,.0);
	complex <float> cn_lastbin(.0,.0);
	complex <float> cdt(.0,.0);
	complex <float> cdw(.0,.0);
	complex <float> norm(1.0,.0);
	complex <float> one(1.0,.0);
	complex <float> two(2.0,.0);
	/*END COMPLEX <FLOAT> DECLARATIONS*/

	/*ARRAY DECLARATIONS*/
	complex <float> *waveletArray;
	complex <float> **datapoint;//original data
	complex <float> **modifieddatapoint;//modification to data e.g. mean-centred data
	spectrum **FT;//FT spectrum
	complex <float> **WT;//WT spectrum
	complex <float> **noiseConstruct;
	complex <float> **reconstructed;
	/*END ARRAY DECLARATIONS*/

	string patientIndex;
	ostringstream convert4;
	convert4 << patientArray[z];
	patientIndex = convert4.str();

	string pcs;
	cout << "PCs as string:" << endl;
	cin.ignore();
	getline(cin,pcs);

	string PCS;
	ostringstream convert10;
	convert10 << pcs;
	PCS = convert10.str();

	string sourcefilename;
	ostringstream convert5;
	convert5 << "ECGs/P" + patientIndex + "_PC.txt";
	//convert5 << "denoised_" + patientIndex + ".txt";
	//convert5 << "ECGs/Patient " + patientIndex + ", 0-5 seconds, data matrix.txt";
	//convert5 << patientIndex + "/P" + patientIndex + "_R_" + PCS + ".txt";
	sourcefilename = convert5.str();

	string complexTransform;
	ostringstream convert6;
	convert6 << "ECGspectra/complexTransform_" + patientIndex + ".txt";
	complexTransform = convert6.str();

	string FTdenoised;
	ostringstream convert7;
	convert7 << "denoisedECG/denoised_" + patientIndex + ".txt";
	FTdenoised = convert7.str();

	string Tplot;
	ostringstream convert8;
	convert8 << "ECGspectra/plot_" + patientIndex + ".dat";
	Tplot = convert8.str();

	string spectralSpread;
	ostringstream convert9;
	convert9 << "ECGspectra/spectralSpread_leads";
	spectralSpread = convert9.str();

	/*OPEN OFSTREAMS*/
	ofstream transform;
	transform.open(complexTransform.c_str());
	if (transform.fail())
	{
		cerr<<"ERROR opening transform.txt."<<endl;
		exit(1);
	}
	else{}

	ofstream plot;
	plot.open(Tplot.c_str());
	if (plot.fail())
	{
		cerr<<"ERROR opening plot.plt."<<endl;
		exit(1);
	}
	else{}

	ofstream denoised;
	denoised.open(FTdenoised.c_str());
	if (denoised.fail())
	{
		cerr<<"ERROR opening denoised.txt."<<endl;
		exit(1);
	}
	else{}

	ofstream Sspread;
	Sspread.open(spectralSpread.c_str());
	if (Sspread.fail())
	{
		cerr<<"ERROR opening Sspread.txt."<<endl;
		exit(1);
	}
	else{}
	/*END OPEN OFSTREAMS*/

	cout << "Enter source-data filename:" << endl;
	//cin >> sourceDataFilename;

	/*IFSTREAM COUNT DATAPOINTS*/
	ifstream lead;
	lead.open(sourcefilename.c_str());
	if (lead.fail())
	{
		cerr<<"ERROR opening source file."<<endl;
		cout << sourcefilename.c_str() << endl;
		exit(1);
	}
	else{}

	while(!lead.eof())
		{
			getline(lead,input);
			sampleCount = sampleCount + 1;
		}
	sampleCount = sampleCount - 1;//NB THIS VARIES FROM COMPILER TO COMPILER
	lead.close();
	cout << "Number of datapoints in source file: " << sampleCount << endl;
	zeroPads = sampleCount/2;
	cout << "Adding " << zeroPads << "zeropads to ends of ECG." << endl;
	count = 2*zeroPads + sampleCount;
	cout << "Total number of datapoints: " << count << endl;

	/*END IFSTREAM COUNT DATAPOINTS*/

	datapoint = new complex <float> * [leads];
	for(int i = 0; i < leads; i = i + 1)
	{
		datapoint[i] = new complex <float> [count]();
	}

	/*IFSTREAM GET ECG*/
	lead.open(sourcefilename.c_str());
	for(int j(0); j < zeroPads; j = j + 1)
	{
		for(int i = 0; i < leads; i = i + 1)
		{
			datapoint[i][j].real() = .0;
		}
	}
	for(int j(zeroPads);j < zeroPads + sampleCount; j = j + 1)
	{
		getline(lead,input);
		stringstream ss2;
		ss2 << input;
		for(int i = 0; i < leads; i = i + 1)
		{
			ss2 >> datapoint[i][j].real();
		}
	}
	for(int j(zeroPads + sampleCount); j < count; j = j + 1)
	{
		for(int i = 0; i < leads; i = i + 1)
		{
			datapoint[i][j].real() = .0;
		}
	}
	ccount.real() = (float)count;
	lead.close();
	/*END IFSTREAM GET ECG*/

	/*DATA: USER INPUT*/
	cout << "The time interval between data points, dt, is " << dt << " (s)." << endl;
	dt = 0.001;
	cdt.real() = dt;

	n = 1;
	cn.real() = (float)n;

	nyquist = 2/(dt);
	cout << "The Nyquist frequency is " << nyquist << " Hz." << endl;
	cout << "The start frequency, f_i, is 0 (Hz):" << endl;
	f_i = 1;
	w_i = f_i*2*pi.real();

	cout << "The end end frequency, f_f, is " << nyquist << " (Hz)." << endl;
	f_f = 200;
	w_f = f_f*2*pi.real();

	cout << "The frequency resolution, df, is " << 1/(count*dt) << " (Hz)." << endl;

	df = 1.0/(10*count*dt);

	dw = df*2*pi.real();
	cdw.real() = dw;

	freqbins = (int)((w_f-w_i)/dw);

	cout << "df: " << df << endl;
	cout << "df*freqbins: " << df*freqbins << endl;
	/*END DATA: USER INPUT*/

	chooseTransform(&transformChoice);

	if(transformChoice == "W" || transformChoice == "w")
	{
		/*WAVELET PARAMETRISATION*/
		chooseWavelet(&waveletChoice, sds);
		cout << "Segment size: " << n*dt << " s." << endl;
		cout << "Enter sigma (s):" << endl;
		//cin >> sigma;
		sigma = 0.05;
		p = (sds*sigma)/(n*dt);
		p = floorf(p + 0.5);
		if (p/2 - (int)(p/2) == 0)
		{
			p = p + 1.0;
			cout << "Modifying number of segments to be odd.  Number of segments: " << (int)p << endl;
		}
		else{}
		cout << "Standard deviation: " << sigma << " s." << endl;
		mu = floorf(p/2.0)*dt;
		cout << "Mean: " << mu << " s." << endl;
		cout << "Building wavelet array." << endl;
		waveletArray = new complex <float> [(int)p]();
		for (int i(0); i < (int)p; i = i + 1)
		{
			if(waveletChoice == "G" || waveletChoice == "g")
			{waveletArray[i].real() = gauss(pi,dt*(float)i,mu,sigma);}
			else if(waveletChoice == "R" || waveletChoice == "r")
			{waveletArray[i].real() = ricker(pi,dt*(float)i,mu,sigma);}
			else if(waveletChoice == "L" || waveletChoice == "l")
			{waveletArray[i].real() = lorentz(pi,dt*(float)i,mu,sigma);cout<<waveletArray[i].real()<<endl;}
			else if(waveletChoice == "H" || waveletChoice == "h")
			{waveletArray[i].real() = haar(i,p);}
			else if(waveletChoice == "K" || waveletChoice == "k")
			{waveletArray[i].real() = 1;}
			else
			{
				cerr<<"ERROR selecting wavelet."<<endl;
				exit(1);
			}
		}
		//plot << "#time\t";
		for(int i(0); i < leads-1; i = i + 1)
		{
			//plot << "lead" << i + 1 << "\t";
		}
		//plot << "lead" << leads << endl;
		/*END WAVELET PARAMETRISATION*/
	}
	else if(transformChoice == "F" || transformChoice == "f")
	{
		/*FOURIER PARAMETRISATION*/
		plot << "#f\t";
		for(int i(0); i < leads-1; i = i + 1)
		{
			plot << "lead" << i + 1 << "\t";
		}
		plot << "lead" << leads << endl;
		/*END FOURIER PARAMETRISATION*/
	}
	else
	{
		cerr<<"ERROR selecting transform choice."<<endl;
		exit(1);
	}

	/*DECLARE ARRAYS*/
	determineBins(count, n, &timebins, &n_lastbin);
	cbins.real() = float(timebins);
	cn_lastbin.real() = float(n_lastbin);

	FT = new spectrum * [leads];
	for(int i = 0; i < leads; i = i + 1)
	{
		FT[i] = new spectrum [freqbins]();
	}

	modifieddatapoint = new complex <float> * [leads];
	WT = new complex <float> * [leads];
	for(int i = 0; i < leads; i = i + 1)
	{
		modifieddatapoint[i] = new complex <float> [timebins]();
		WT[i] = new complex <float> [timebins]();
	}
	/*END DECLARE ARRAYS*/
	meanCentre(datapoint,count,leads);
	/*LEAD LOOP*/
	for(int leadNo(0); leadNo < leads; leadNo = leadNo + 1)
	{
		if(n==1)
		{
				for(int j(0); j < count; j = j + 1)
				{modifieddatapoint[leadNo][j] = datapoint[leadNo][j];}
		}
		else //calculates mean functions each over n points, and assigns to new modifieddatapoint with 'timebins' elements
		{
			for(int i(0);i<(timebins-1)*n;i=i+n)
			{
				for(int j=i;j<i+n;j=j+1)
				{
					sum = sum + datapoint[leadNo][j];
				}
				modifieddatapoint[leadNo][(i/n)] = sum/cn;
				sum = (.0,.0);
			}

			for(int i=(timebins-1)*n;i<count;i=i+1)
			{
				sum = sum + datapoint[leadNo][i];
			}
			modifieddatapoint[leadNo][timebins-1] = sum/(ccount-(cbins-one)*cn);
			sum = (.0,.0);
		}

		if(transformChoice == "W" || transformChoice == "w")
		{
		//for() FOR LOOP TO MODIFY PARAMETERS OF WAVELET
		//{
			int mMin(0);
			int mMax(0);
			int nMin(0);
			int nMax(0);
			float sumWT(0.0);
			for(int k(0);k<timebins;k=k+1) //FOR LOOP TO MODIFY POSITION OF WAVELET RELATIVE TO LEAD
			{
				waveletRange(p,k,timebins,&mMin,&mMax);
				for(int m(mMin); m <= mMax; m = m + 1) //FOR LOOP TO CALCULATE CONVOLUTION AT EACH POSITION
				{
					clocalFunction = modifieddatapoint[leadNo][k-(int)floorf(p/2.0)+m];
					WT[leadNo][k] = WT[leadNo][k] + clocalFunction*waveletArray[m]*cn*cdt/norm;
				}
				sumWT = sumWT + WT[leadNo][k].real();
			}
			cout << sumWT << endl;
			sumWT = 0.0;
		//}
		}
		else if(transformChoice == "F" || transformChoice == "f")
		{
			for(int freqbin(0); freqbin < freqbins; freqbin = freqbin + 1)
			{
				for(int timebin(0); timebin < timebins; timebin = timebin + 1)
				{
					arg = -(freqbin + 1)*dw*(timebin + 1)*dt;
					expo.real() = cos(arg);
					expo.imag() = sin(arg);
					FT[leadNo][freqbin] = FT[leadNo][freqbin] + modifieddatapoint[leadNo][timebin]*expo;
				}
			}
		}
		else
		{
			cerr << "ERROR selecting transform choice." << endl;
			exit(1);
		}
		cout << "Patient " << patientIndex << ", Lead " << leadNo + 1 << " transform completed." << endl;

	}
/*END LEAD LOOP*/

	/*PRINT DATA*/
	if(transformChoice == "W" || transformChoice == "w")
	{
		for(int i(0);i < timebins;i = i + 1)
		{
			transform << i*n*dt << "\t";
			//plot << i*n*dt << "\t";
			for(int leadNo(0); leadNo < leads; leadNo = leadNo + 1)
			{
				if(leadNo < leads)
				{
					transform << WT[leadNo][i] << "\t";
					plot << (WT[leadNo][i]).real() << "\t";//************************************************************
				}
				else if(leadNo == leads-1)
				{
					transform << WT[leadNo][i];
					plot << (WT[leadNo][i]).real();//************************************************************
				}
				else
				{
					cerr << "ERROR: leadNo > leads." << endl;
					exit(1);
				}
			}
			transform << endl;
			plot << endl;
			
		}
		cout << "Patient " << patientIndex << ". Print completed." << endl;
	}
	else if(transformChoice == "F" || transformChoice == "f")
	{
		for(int leadNo(0); leadNo < leads; leadNo = leadNo + 1)
		{
			for(int freqbin(0); freqbin < freqbins; freqbin = freqbin + 1)
			{
				sumFT = sumFT + FT[leadNo][freqbin].a2();
				sumFT2 = sumFT2 + pow(FT[leadNo][freqbin].a2(),2.0);
			}
			spread = pow((sumFT2/freqbins - pow(sumFT/freqbins,2.0)),-0.5);
			sumFT2 = sumFT = 0.0;
			//cout << spread << endl;
		}

		for(int freqbin(0); freqbin < freqbins; freqbin = freqbin + 1)
		{
			transform << (freqbin+1)*df << "\t";
			plot << (freqbin+1)*df << "\t";
			for(int leadNo(0); leadNo < leads; leadNo = leadNo + 1)
			{
				if(leadNo < leads - 1)
				{
					transform << FT[leadNo][freqbin].FT() << "\t";
					plot << FT[leadNo][freqbin].a2() << "\t";
				}
				else if(leadNo == leads-1)
				{
					transform << FT[leadNo][freqbin].FT();
					plot << FT[leadNo][freqbin].a2();
				}
				else
				{
					cerr << "ERROR: leadNo > leads." << endl;
					exit(1);
				}
			}
			transform << endl;
			plot << endl;
		}
		cout << "Patient " << patientIndex << ". Print completed." << endl; 
	}
	else{}
	/*END PRINT DATA*/

	/*CLOSE*/
	plot.close();
	transform.close();
	/*END CLOSE FILES*/

	/*PLOT DATA*/
	cout << "[P]lot leads using gnuplot?" << endl;
	//cin >> plotChoice;
	plotChoice = "P";
	if(plotChoice == "P" || plotChoice == "p")
	{
		if(transformChoice == "W" || transformChoice == "w")
		{
			writeWTplt(leads,patientIndex);
			//system("gnuplot WT.plt");
		}
		else if(transformChoice == "F" || transformChoice == "f")
		{
			writeFTplt(leads,patientIndex);
			system("gnuplot FT.plt");
		}
		else
		{
			cerr << "ERROR selecting transform choice." << endl;
			exit(1);
		}
	}
	else{}
	/*END PLOT DATA*/

	/*THRESHOLDING*/
	cout << "[D]enoise leads?" << endl;
	//cin >> denoiseChoice;
	denoiseChoice = "n";
	if(transformChoice == "F" || transformChoice == "f")
	{
	if(denoiseChoice == "D" || denoiseChoice == "d")
	{
		noiseConstruct = new complex <float> * [leads];
		reconstructed = new complex <float> * [leads];
		for(int i = 0; i < leads; i = i + 1)
		{
			noiseConstruct[i] = new complex <float> [timebins]();
			reconstructed[i] = new complex <float> [timebins]();
		}
		//write function that takes start and end frequency as parameters
		float lBoundW = 2*pi.real()*49;
		float uBoundW = 2*pi.real()*500;
		//float lBoundW = 310.0;
		//float uBoundW = 320.0;
		for(int leadNo(0); leadNo < leads; leadNo = leadNo + 1)
		{
			for(int freqbin(0); freqbin < freqbins; freqbin = freqbin + 1)
			{
				if((freqbin + 1)*dw >= lBoundW && (freqbin + 1)*dw <= uBoundW)
				{
					for(int timebin(0); timebin < timebins; timebin = timebin + 1)
					{
						
						arg = timebin*dt*(freqbin + 1)/freqbins;
						expo.real() = cos(arg);
						expo.imag() = sin(arg);
						noiseConstruct[leadNo][timebin] = noiseConstruct[leadNo][timebin] + FT[leadNo][freqbin].FT()*expo/sqrt((float)freqbins);
					}
				}
				else if((freqbin + 1)*dw <= lBoundW || (freqbin + 1)*dw >= uBoundW)
				{
				}
				else
				{
					cerr << "ERROR: freqbin out of bounds." << endl;
					exit(1);
				}
			}
		}
		
		for(int leadNo(0); leadNo < leads; leadNo = leadNo + 1)
		{
			
			for(int k(0); k < timebins; k = k + 1)
			{
				reconstructed[leadNo][k] = datapoint[leadNo][k] - noiseConstruct[leadNo][k];
			}
		}
		for(int i(0);i < timebins;i = i + 1)
		{
			//denoised << i*n*dt << "\t";
			for(int leadNo(0); leadNo < leads; leadNo = leadNo + 1)
			{
				if(leadNo < leads)
				{
					denoised << reconstructed[leadNo][i].real() << "\t";//*****************************************************
				}
				else if(leadNo == leads-1)
				{
					denoised << reconstructed[leadNo][i].real();//************************************************************
				}
				else
				{
					cerr << "ERROR: leadNo > leads." << endl;
					exit(1);}
				}
			denoised << endl;
			
		}
		cout << "Patient " << patientIndex << ". Denoised signal written to denoisedECG/denoised_" + patientIndex + ".txt";

	/*DELETE*/
	denoised.close();
	/*2D ARRAYS*/

	for(int i(0); i < leads; i = i + 1)
	{
		string command1;
		string command2;
		ostringstream convert9;
		ostringstream convert10;
		convert10 << i + 1;
		command1 = convert10.str();
		//convert9 << "convert ECGspectra/FT_lead" + command1 + "_patient" + patientIndex + ".svg ECGspectra/FT_lead_" + command1 + "_patient" + patientIndex + ".jpg";
		convert9 << "convert ECGspectra/FT_lead" + command1 + "_patient" + patientIndex + ".svg ECGspectra/" + command1 + ".jpg";
		command2 = convert9.str();
		system(command2.c_str());
	}
	cout << "Conversion to .jpg complete." << endl;
	}
	else{}
	}
	else if(transformChoice == "W" || transformChoice == "w")
	{
		/*for(int i(0); i < leads; i = i + 1)
	{
		string command1;
		string command2;
		ostringstream convert9;
		ostringstream convert10;
		convert10 << i + 1;
		command1 = convert10.str();
		convert9 << "convert ECGspectra/WT_lead" + command1 + "_patient" + patientIndex + ".svg ECGspectra/FT_lead_" + command1 + "_patient" + patientIndex + ".jpg";
		command2 = convert9.str();
		system(command2.c_str());
	}
	cout << "Conversion to .jpg complete." << endl;*/
	}
	else
	{
		cerr << "ERROR selecting transform." << endl;
		exit(1);
	}
	/*END THRESHOLDING*/

	//spreadArray = new float * [patients];
	for(int patient = 0; patient < patients; patient = patient + 1)
	{
		delete [] spreadArray[patient];
	}
	delete [] spreadArray;

	for(int i = 0; i < leads; i = i + 1)
	{
		delete [] datapoint[i];
		delete [] modifieddatapoint[i];
		delete [] WT[i];
		delete [] FT[i];
		if(transformChoice == "D" || transformChoice == "d")
		{
			delete [] noiseConstruct;
			delete [] reconstructed;
		}
	}
	delete [] datapoint;
	delete [] modifieddatapoint;
	delete [] FT;
	delete [] WT;
	if(denoiseChoice == "D" || denoiseChoice == "d")
	{
		delete [] noiseConstruct;
		delete [] reconstructed;
	}
	/*END 2D ARRAYS*/

	/*1D ARRAYS*/
	if(transformChoice == "W" || transformChoice == "w")
	{delete [] waveletArray;}
	else if(transformChoice == "F" || transformChoice == "f"){}
	else
	{
		cerr << "ERROR selecting transform." << endl;
		exit(1);
	}
	/*END 1D ARRAYS*/

	/*END DELETE*/
}
/*END PATIENT LOOP*/
	delete [] patientArray;
	cout << "Dynamic arrays deleted." << endl;
	return(0);
}

void t(int i)
{cout<<"test "<<i<<endl;}

/*MEAN*/
void meanCentre(complex <float> ** datapoint, const int count, const int leads)
{
	for(int leadNo(0); leadNo < leads; leadNo = leadNo + 1)//sets baseline to mean of leadArray i.e. new baseline at zero
	{
		complex <float> sum(.0,.0);
		complex <float> mean(.0,.0);
		for(int j(0); j < count; j = j + 1)
		{sum = sum + datapoint[leadNo][j];}
		mean = sum/(complex <float>)count;
		for(int j(0); j < count; j = j + 1)
		{datapoint[leadNo][j] = datapoint[leadNo][j] - mean;}
		sum = (.0,.0);
		mean = (.0,.0);
	}
}
/*END MEAN*/

void writeWTplt(int leads,string patientIndex)
{
	ofstream myFile1;
	myFile1.open("WT.plt");
	myFile1 << "#reset #reset gnuplot options" << endl;
	myFile1 << "set terminal svg enhanced size 4000 3000 fname \"Times\" fsize 36" << endl;
	for(int i(0); i < leads; i = i + 1)
	{
		string colNo;
		string fileNo;
		ostringstream convert1;
		ostringstream convert2;
		convert1 << (i + 2);
		convert2 << (i + 1);
		colNo = convert1.str();
		fileNo = convert2.str();
		myFile1 << "set output \"ECGspectra/WT_lead" + fileNo + "_patient" + patientIndex + ".svg\"" << endl;
		myFile1 << "set title \"Wavelet, lead " + fileNo + ", patient " + patientIndex + "\"" << endl;
		myFile1 << "set xlabel \"Time (s)\"" << endl;
		myFile1 << "set ylabel \"Convolution (mV s)\"" << endl;
		myFile1 << "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb\"#fbffc1\" behind" << endl;
		//myFile1 << "plot x linecolor rgb \"blue\" linewidth 2 pointtype 4 pointsize 2" << endl;
		myFile1 << "plot \"./ECGspectra/plot_" + patientIndex + ".dat\" using 1:" + colNo + " title \"\" pointtype 7 pointsize 1" << endl;
	}
	myFile1.close();
}

void writeFTplt(int leads, string patientIndex)
{
	ofstream myFile1;
	myFile1.open("FT.plt");
	myFile1 << "#reset #reset gnuplot options" << endl;
	myFile1 << "set terminal svg enhanced size 4000 3000 fname \"Times\" fsize 36" << endl;
	for(int i(0); i < leads; i = i + 1)
	{
		string colNo;
		string fileNo;
		ostringstream convert1;
		ostringstream convert2;
		convert1 << (i + 2);
		convert2 << (i + 1);
		colNo = convert1.str();
		fileNo = convert2.str();
		myFile1 << "set output \"ECGspectra/FT_lead" + fileNo + "_patient" + patientIndex + ".svg\"" << endl;
		myFile1 << "set title \"Fourier, lead " + fileNo + ", patient " + patientIndex + "\"" << endl;
		myFile1 << "set xlabel \"Frequency (Hz)\"" << endl;
		myFile1 << "set ylabel \"Convolution (mV s)\"" << endl;
		myFile1 << "set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb\"#fbffc1\" behind" << endl;
		//myFile1 << "plot x linecolor rgb \"blue\" linewidth 2 pointtype 4 pointsize 2" << endl;
		myFile1 << "plot \"./ECGspectra/plot_" + patientIndex + ".dat\" using 1:" + colNo + " title \"\" pointtype 7 pointsize 1" << endl;
	}
	myFile1.close();
}

void waveletRange(float p, int k, int timebins, int *mMin, int *mMax)
{
	if((int)floorf(p/2.0) - k >= 0)
	{
		*mMin = (int)floorf(p/2.0) - k;
	}
	else
	{
		*mMin = 0;
	}
		if(timebins - k + (int)floorf(p/2.0) <= p - 1)
	{
		*mMax = timebins - k + (int)floorf(p/2.0) - 1;
	}
	else
	{
		*mMax = p - 1;
	}
}

void chooseWavelet(string *waveletChoice, const float sds)
{
	bool binary(1);
	do
	{
		cout << "Select wavelet: [G]aussian, [R]icker [L]orentzian, [H]aar, [K]ronecker-Delta." << endl;
		//cin >> *waveletChoice;
		*waveletChoice = "G";
		if(*waveletChoice == "G" || *waveletChoice == "g")
		{
			binary = 1;
			cout << "Using Gaussian wavelet function, psi(t,mu,sigma) = A*exp(-(t-mu)^2/(2*sigma^2)), on the range [-"<<sds<<" sigma,+"<<sds<<" sigma]." << endl;
		}
		else if(*waveletChoice == "R" || *waveletChoice == "r")
		{
			binary = 1;
			cout << "Using the Ricker wavelet function, psi(t) = 2/(3^0.5*pi^0.25)*(1-(t-mu)^2/sigma^2)*exp(-(t-mu)^2/(2*sigma^2)), on the range [-"<<sds<<" sigma,+"<<sds<<" sigma]." << endl;
		}
		else if(*waveletChoice == "L" || *waveletChoice == "l")
		{
			binary = 1;
			cout << "Using the Lorentzian wavelet function, psi(t) = ***, on the range [-"<<sds<<" sigma,+"<<sds<<" sigma]." << endl;
		}
		else if(*waveletChoice == "H" || *waveletChoice == "h")
		{
			binary = 1;
			cout << "Using the Haar wavelet function, psi(t) = ***, on the range [-"<<sds<<" sigma,+"<<sds<<" sigma]." << endl;
		}
		else if(*waveletChoice == "K" || *waveletChoice == "k")
		{
			binary = 1;
			cout << "Using the Kronecker-Delta wavelet function, psi(t) = ***, on the range [-"<<sds<<" sigma,+"<<sds<<" sigma]." << endl;
		}
		else
		{
			binary = 0;
			cout << "Not a valid wavelet." << endl;
		}
	}
	while(!binary);
}

void determineBins(int count, int n, int * timebins, int * n_lastbin)
{
	if((count/n)*n<count)
	{
		*timebins = count/n + 1;
		*n_lastbin = count - (count/n)*n;
		cout << "The last bin contains " << *n_lastbin << " datapoints." << endl;
	}
	else if((count/n)*n==count)
	{
		*timebins = count/n;
		cout << "The last bin is full." << endl;
	}
	else
	{
		cout << "ERROR: bins." << endl;
		exit(1);
	}	
	cout << "There are " << *timebins << " time bins." << endl;
}

void chooseTransform(string *transformChoice)
{
	bool binary(1);
	do
	{
		cout << "Do you wish to perform a [F]ourier or [W]avelet transform?" << endl;
		//cin >> *transformChoice;
		*transformChoice = "F";
		if(*transformChoice == "F" || *transformChoice == "f")
		{
			binary = 1;
			cout << "Using Fourier." << endl;
		}
		else if(*transformChoice == "W" || *transformChoice == "w")
		{
			binary = 1;
			cout << "Using Wavelet." << endl;
		}
		else
		{
			binary = 0;
			cout << "Not a valid transform." << endl;
		}
	}
	while(!binary);
}

/*WAVELET DEFINITIONS*/
float gauss(const complex <float> pi, const float x, const float mu, const float sigma)
{return exp(-pow((x-mu),2.0)/(2.0*pow(sigma,2.0)));}
float ricker(const complex <float> pi, const float x, const float mu, const float sigma)
{return 2.0/(pow(3.0,0.5)*pow(pi.real(),0.25))*(1.0-pow(x-mu,2.0)/pow(sigma,2.0))*exp(-pow((x-mu),2.0)/(2.0*pow(sigma,2.0)));}
float lorentz(const complex <float> pi, const float x, const float mu, const float sigma)
{return 15.0*sigma*pow(pi.real(),-1.0)*pow(pow((x-mu),2.0)+pow(0.5*sigma,2.0),-1.0);}
float haar(const int i, const float p)
{
	if(i<(int)floorf(p/2.0))
		{return 1.0;}
		else if(i==(int)floorf(p/2.0))
		{return 0.0;}
		else if(i>(int)floorf(p/2.0))
		{return -1.0;}
}
/*END WAVELET DEFINITIONS*/
