#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <utility>
#include <cassert>
#include "ec_ops.h"
using namespace std;

Zp Zp::inverse() const{
	// Implement the Extended Euclidean Algorithm to return the inverse mod PRIME		
	uberzahl u1(PRIME);
	uberzahl u2(value);
	uberzahl q, t1(0), t2(1), temp;
	while(u2 != 0){
		temp = u2;
		q = u1 / u2;
		u2 = u1 - q * u2;
		u1 = temp;

		temp = t2;
		t2 = t1 - q * t2;
		t1 = temp;
	}
	return Zp(t1);
}

ECpoint ECpoint::operator + (const ECpoint &a) const {
	// Implement  elliptic curve addition 		
	// If either one is infinity point
	if(a.infinityPoint){
		return *this;
	}
	else if(infinityPoint){
		return a;
	}

	Zp slope, rx, ry;
	Zp px = x, py = y, 
	   qx = a.x, qy = a.y;
	// If P != Q, PX != QX
	if(!(px == qx)){
		slope = (qy - py) * (qx - px).inverse();
	}
	// If P = Q, py != 0
	else if(!(py == Zp(0))){
		slope = ((px * px * Zp(3)) + Zp(A)) * (py + py).inverse();
	}
	// Other cases, P + Q = O
	else{
		return ECpoint(true);
	}
	rx = slope * slope - px - qx;
	ry = slope * (px - rx) - py;
	return ECpoint(rx, ry); 
}

ECpoint ECpoint::repeatSum(ECpoint p, uberzahl v) const {
	//Find the sum of p+p+...+p (vtimes)
	// If P = O
	if(p.infinityPoint){
		return p;
	}
	// If v >= EC order
	while(v >= ORDER){
		v = v - ORDER;	
	}
	if(v == 0){
		return ECpoint(true);
	}
	
	ECpoint f(true);
	for(unsigned expOrder = v.bitLength(); expOrder > 0; --expOrder){
		f = f + f;
		if(v.bit(expOrder - 1) == 1){
			f = f + p; 
		}
	}
	return f;
}

Zp ECsystem::power(Zp val, uberzahl pow) {
	//Find the product of val*val*...*val (pow times)
	// If val = 0 or val = 1
	if(val == Zp(0) || val == Zp(1)){
		return val;
	}
	while(pow >= PRIME){
		pow = pow - PRIME;
	}
	// If pow = 0;
	if(pow == 0){
		return Zp(1);
	}

	Zp f(1);
	for(unsigned expOrder = pow.bitLength(); expOrder > 0; --expOrder){
		f = f * f;
		if(pow.bit(expOrder - 1) == 1){
			f = f * val; 
		}
	}
	return f;
}


uberzahl ECsystem::pointCompress(ECpoint e) {
	//It is the gamma function explained in the assignment.
	//Note: Here return type is mpz_class because the function may
	//map to a value greater than the defined PRIME number (i.e, range of Zp)
	//This function is fully defined.	
	uberzahl compressedPoint = e.x.getValue();
	compressedPoint = compressedPoint<<1;
	
	if(e.infinityPoint) {
		cout<<"Point cannot be compressed as its INF-POINT"<<flush;
		abort();
		}
	else {
		if (e.y.getValue().bit(0) == 1)
			compressedPoint = compressedPoint + 1;
		}
	//cout<<"For point  "<<e<<"  Compressed point is <<"<<compressedPoint<<"\n";
	return compressedPoint;

}

ECpoint ECsystem::pointDecompress(uberzahl compressedPoint){
	//Implement the delta function for decompressing the compressed point
 	Zp rx(compressedPoint >> 1), ry, z;
	z = (rx * rx * rx) + (Zp(A) * rx) + Zp(B);
	ry = power(z, (PRIME + 1) >> 2);
	if(ry.getValue().bit(0) == compressedPoint.bit(0)){
		return ECpoint(rx, ry);
	}
	else{
		return ECpoint(rx, Zp(PRIME) - ry);
	}
}


pair<pair<Zp,Zp>,uberzahl> ECsystem::encrypt(ECpoint publicKey, uberzahl privateKey, Zp plaintext0, Zp plaintext1){
	// You must implement elliptic curve encryption
	// Do not generate a random key. Use the private key that is passed from the main function
	// 1. Generate random XB
	// 2. Compute Q = XB * G; R = XB * P = XB * XA * G
	ECpoint Q = privateKey * G, 
			R = privateKey * publicKey;
	// 3. Set C0 = M0 * Rx; C1 = M1 * Ry; C2 = gamma(Q)
	Zp C0 = Zp(plaintext0) * R.x, 
	   C1 = Zp(plaintext1) * R.y;
	uberzahl C2 = pointCompress(Q);
	return make_pair(make_pair(C0, C1), C2);
}


pair<Zp,Zp> ECsystem::decrypt(pair<pair<Zp,Zp>, uberzahl> ciphertext){
	// Implement EC Decryption
	ECpoint R = privateKey * pointDecompress(ciphertext.second);
	pair<Zp, Zp> plaintext;
	plaintext.first = R.x.inverse() * ciphertext.first.first;
	plaintext.second = R.y.inverse() * ciphertext.first.second;
	return plaintext;
}


/*
 * main: Compute a pair of public key and private key
 *       Generate plaintext (m1, m2)
 *       Encrypt plaintext using elliptic curve encryption
 *       Decrypt ciphertext using elliptic curve decryption
 *       Should get the original plaintext
 *       Don't change anything in main.  We will use this to 
 *       evaluate the correctness of your program.
 */


int main(void){
	srand(time(0));
	ECsystem ec;
	uberzahl m0, m1, xa, xb;
	m0 = random("0", PRIME);	
	m1 = random("0", PRIME);	
	xa = random("0", PRIME);	
	xb = random("0", PRIME);	
	pair <ECpoint, uberzahl> keys = ec.generateKeys(xa);
	Zp plaintext0(m0);
	Zp plaintext1(m1);
	ECpoint publicKey = keys.first;
	pair<pair<Zp,Zp>, uberzahl> ciphertext = ec.encrypt(publicKey, xb, plaintext0, plaintext1);	
	pair<Zp,Zp> plaintext_out = ec.decrypt(ciphertext);	
	if(plaintext0 == plaintext_out.first && plaintext1 == plaintext_out.second){
		cout << "Correct!" << endl;
	}
	else{
		cout << "Private key(Sender, Receiver): (" << xb << ", " << xa << "); Public key: " << publicKey << endl;
		cout << "Original plaintext: (" << plaintext0 << ", " << plaintext1 << ")\n";
		cout << "Decrypted plaintext: (" << plaintext_out.first << ", " << plaintext_out.second << ")\n";
	}
	return 1;
}


