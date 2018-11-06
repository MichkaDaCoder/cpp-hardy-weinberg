#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace std;

/**
* HardyWeinberg.cpp :
* @author: MichkaDaCoder
* @since : 27.10.2018
*
* Compute Alleles frequencies, homozygote/heterozygotes frequencies,
* expected genetic number of alleles and Chi Square value using Hardy Weinberg law.
*
* Genetic is awesome. :)
* Thank you.
**/


class HardyWeinberg {
private:
  // AA genotype
  double AA;
  // Aa genotype
  double Aa;
  // aa genotype
  double aa;

public :

    ///Constructor Logic
    HardyWeinberg(double AA,double Aa, double aa) {
       this->AA=AA;
       this->Aa=Aa;
       this->aa=aa;
    }

    ~HardyWeinberg() {}

    ///Getters
    double getAA() {return this->AA;}
    double getAa() {return this->Aa;}
    double getaa() {return this->aa;}

    ///Setters
    void setAA(double AA) {this->AA=AA;}
    void setAa(double Aa) {this->Aa=Aa;}
    void setaa(double aa) {this->aa=aa;}

    /** getSum
        Compute the sum of the Alleles numbers
        @param AA
        @param Aa
        @param aa
        @return double sum of the allele numbers.
        @author MichkaDaCoder
    */
    double getSum(double AA, double Aa,double aa) {
       return AA+Aa+aa;
    }

    /** getP
        Compute the frequency of the dominant Allele 'A'
        @param AA
        @param Aa
        @param aa
        @return double frequency of the dominant Allele 'A'
        @author MichkaDaCoder
     */
    double getP(double AA,double Aa, double aa) {
      return ((2*AA)+Aa)/(2*(AA+Aa+aa));
    }

    /** getQ
        Compute the frequency of the recessive Allele 'a'
        @param p
        @return double frequency of the recessive Allele 'a'
        @author MichkaDaCoder
     */
    double getQ(double p) {
      return 1-p;
    }

    /** getAAFrequency
        Compute the frequence of the homozygote genotype 'AA'
        @param p
        @return double frequence of AA
        @author MichkaDaCoder
     */
    double getAAFrequency(double p) {
      return p*p;
    }

     /** getAaFrequency
        Compute the frequence of the heterozygote genotype 'Aa'
        @param p
        @param q
        @return double frequence of Aa
        @author MichkaDaCoder
     */
    double getAaFrequency(double p, double q) {
      return 2*p*q;
    }

    /** getaaFrequency
        Compute the frequency of the homozygote genotype 'aa'
        @param q
        @return number double
        @author MichkaDaCoder
    */
    double getaaFrequency(double q) {
      return q*q;
    }

    /**getExpAA
       Compute the expected genetype number of the homozygote allele 'AA'
       @param p
       @param sum
       @return sum double
       @author MichkaDaCoder
    */
    double getExpAA(double p, double sum) {
       return (p*p)*sum;
    }

    /**getExpAa
       Compute the expected genetype number of the heterogygote allele Aa
       @param p
       @param q
       @param sum
       @return number double
       @author MichkaDaCoder
    */
    double getExpAa(double p, double q, double sum) {
       return 2*p*q*sum;
    }

    /**getExpaa
       Compute the expected genetype number of the homozygote allele 'aa'
       @param q
       @param sum
       @return number double
       @author MichkaDaCoder
    */
    double getExpaa(double q, double sum) {
       return (q*q)*sum;
    }

    /**getKhiSquare
       Compute the expected genetype number of the homozygote allele 'aa'
       @param q
       @param sum
       @return number double
       @author MichkaDaCoder
    */
    double getKhiSquare(double AA, double expAA, double Aa, double expAa, double aa, double expaa) {
       return (pow((AA-expAA),2)/(expAA)) + (pow((Aa-expAa),2)/(expAa)) + (pow((aa-expaa),2)/(expaa));
    }
};

/* Entry point */
int main() {
system("title Hardy Weinberg Principal"); //!Setting a title to the console
double AA,Aa,aa;
cout<<"Enter number of Phenotype AA: ";
cin>>AA;
cout<<"Enter number of Phenotype Aa: ";
cin>>Aa;
cout<<"Enter number of Phenotype aa: ";
cin>>aa;
HardyWeinberg*hw=new HardyWeinberg(AA,Aa,aa);
double p=hw->getP(hw->getAA(),hw->getAa(),hw->getaa());
double q=hw->getQ(hw->getP(hw->getAA(),hw->getAa(),hw->getaa()));
double p2=hw->getAAFrequency(hw->getP(hw->getAA(),hw->getAa(),hw->getaa()));
double pq2=hw->getAaFrequency(hw->getP(hw->getAA(),hw->getAa(),hw->getaa()),hw->getQ(hw->getP(hw->getAA(),hw->getAa(),hw->getaa())));
double q2=hw->getAaFrequency(hw->getP(hw->getAA(),hw->getAa(),hw->getaa()),hw->getQ(hw->getP(hw->getAA(),hw->getAa(),hw->getaa())));
double sum=hw->getSum(AA,Aa,aa);
double expAA=hw->getExpAA(p, sum);
double expAa=hw->getExpAa(p,q,sum);
double expaa=hw->getExpaa(q,sum);
double khi=hw->getKhiSquare(AA,expAA,Aa, expAa, aa, expaa);
double fAA=hw->getAAFrequency(hw->getP(hw->getAA(),hw->getAa(),hw->getaa()));
double fAa=hw->getAaFrequency(hw->getP(hw->getAA(),hw->getAa(),hw->getaa()),hw->getQ(hw->getP(hw->getAA(),hw->getAa(),hw->getaa())));
double faa=hw->getQ(p);
system("cls");
cout<<"Phenotype["<<hw->getAA()<<","<<hw->getAa()<<","<<hw->getaa()<<"]"<<endl;
cout<<"Total: "<<sum<<endl<<endl;

cout<<"Computing p and q..."<<endl<<endl;
cout<<"p= "<<p<<endl;
cout<<"q= "<<q<<endl<<endl;

cout<<"Computing Alleles Frequencies..."<<endl<<endl;
cout<<"F(AA)= "<<fAA<<endl;
cout<<"F(Aa)= "<<fAa<<endl;
cout<<"F(a)= "<<faa<<endl;
cout<<"p²+2pq+q²= "<<p2+pq2+q2<<endl<<endl;
cout<<"Computing Expected Genotype Numbers..."<<endl<<endl;
cout<<"Exp(AA)= "<<expAA<<endl;
cout<<"Exp(Aa)= "<<expAa<<endl;
cout<<"Exp(aa)= "<<expaa<<endl<<endl;
cout<<"Computing Khi Square..."<<endl;
cout<<"Khi Square= "<<khi<<endl<<endl;

cin.get();
cin.get();

system("cls");
cout<<" Brought to you by : MichkaDaCoder"<<endl;
cout<<" github.com/michkadacoder"<<endl<<endl;
cout<<" Thank you and good bye :)"<<endl<<endl;
cin.ignore();
cin.get();

}
