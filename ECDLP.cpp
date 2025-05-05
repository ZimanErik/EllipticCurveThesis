#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <iostream>
#include <exception>
#include <sstream>
#include <chrono>
#include <unordered_map>
#include <streambuf>
#include <vector>


#include <fstream>
#include <vector>
#include <string>
#include <nlohmann/json.hpp>

NTL_CLIENT

class NullStreamBuf : public std::streambuf {
    int overflow(int ch) override {
        return ch; 
    }
};
   
std::streambuf* originalCoutBuf = std::cout.rdbuf();

void turnOffCout() {
    static NullStreamBuf nullStreamBuf;  
    std::cout.rdbuf(&nullStreamBuf);
}

void turnOnCout() {
    std::cout.rdbuf(originalCoutBuf);
}
   
bool isPrime(ZZ n) {
    cout << " started checking primality of n :" << n << endl;
    if (n <= ZZ(1)) {
      return false;
    }
    for (ZZ i = ZZ(2); i * i <= n; ++i) {
      cout << "still chekicng";
      if (n % i == 0) {
        return false;
      }
    }
    return true;
}

struct ECPoint {
    ZZ x, y;
    bool isPointAtInfinity;

    ECPoint() : isPointAtInfinity(true) {}  
    ECPoint(ZZ x, ZZ y) : x(x), y(y), isPointAtInfinity(false) {}

    void print() const {
        if (isPointAtInfinity) cout << "IPoint(∞)" << endl;
        else cout << "(" << x << ", " << y << ")" << endl;
    }

    string toString() const {
        ostringstream oss;
        if (isPointAtInfinity)
            oss << "IPoint(∞)";
        else
            oss << "(" << x << ", " << y << ")";
        return oss.str();
    }

    bool operator==(const ECPoint& other) const {
        if (isPointAtInfinity && other.isPointAtInfinity) return true;
        if (isPointAtInfinity || other.isPointAtInfinity) return false;
        return (x == other.x) && (y == other.y);
    }
};

struct ZZHasher {
    size_t operator()(const ZZ& num) const {
        stringstream ss;
        ss << num;  // Convert ZZ to string
        string str = ss.str();
        return std::hash<string>{}(str);  // Hash the string representation of ZZ
    }
};

// (1,23) (12,3) TODO check if these hashes are equal 
struct ECPointHasher
{
  size_t operator()(const ECPoint& P) const
  {
    size_t rowHash = ZZHasher{}(P.x);
    size_t colHash = ZZHasher{}(P.y) << 1;
    return rowHash ^ colHash;
  }
};


struct EllipticCurve {
    ZZ a, b, p;

    EllipticCurve(ZZ a, ZZ b, ZZ p) : a(a), b(b), p(p) {}

    // TODO : check if there is a way to calculate moduli before and restore it after function ends     
    bool isOnCurve(const ECPoint& P) const {

        if (P.isPointAtInfinity) return true;
        ZZ_p::init(p);
        
        ZZ_p x = to_ZZ_p(P.x);
        ZZ_p y = to_ZZ_p(P.y);
        
        ZZ_p left = sqr(y);
        
        //? IDEA:  ZZ_p right = x*x*x + conv<ZZ_p>(a) * x + conv<ZZ_p>(b);
        //? IDEA:  ZZ_p right2 = to_ZZ_p(PowerMod(P.x, to_ZZ(3), p)) + to_ZZ_p(a)*x + to_ZZ_p(b); 
        ZZ_p right = x*x*x + to_ZZ_p(a)*x + to_ZZ_p(b);      
        
        // (y^2) mod p  ==  (x^3 + ax + b) mod p
        // cout << left << "|||" << right << endl;
        return (left == right);
    }

    void print() const {
        cout << "Elliptic curve : y^2 = x^3 + " << a << "x + "<< b << " mod "<< p << endl;
    }

    string toString() const {
        ostringstream oss;
        oss << "Elliptic curve : y^2 = x^3 + " << a << "x + "<< b << " mod "<< p ;
        return oss.str();
    }
};

class PointNotOnCurveException : public exception {
private:
    string message;

public:
    PointNotOnCurveException(const ECPoint& P, const EllipticCurve& EC){
        ostringstream mess;
        mess <<"Point ("<<P.x<<","<<P.y<<") is not on curve y^2=x^3+"<<EC.a<<"x"<<"+"<<EC.b<<"mod "<<EC.p; 
        message = mess.str();
    }

    const char* what() const throw()
    {
        return message.c_str();
    }
};

ECPoint doublePoint(const ECPoint& P,const EllipticCurve& EC){
   
    if(P.y==0) return ECPoint();
    ZZ_p::init(EC.p);
    ZZ_p p_inv = inv(to_ZZ_p(2)*to_ZZ_p(P.y));

    // in reality the second '*' bellow is division because of inverse multiplication 
    ZZ_p m = ( (to_ZZ_p(3) * to_ZZ_p(sqr(P.x)) + to_ZZ_p(EC.a) ) * (p_inv) );    
    ZZ_p n = m * m - to_ZZ_p(2) * to_ZZ_p(P.x);

    ZZ x = rep(n);
    ZZ y = rep( m*(to_ZZ_p(P.x) - n) - to_ZZ_p(P.y));
    return ECPoint(x,y);
}

ECPoint addPoints(const ECPoint& P, const ECPoint& Q, const EllipticCurve& EC){

    try {
        if (!EC.isOnCurve(P)) throw PointNotOnCurveException(P,EC);
        if (!EC.isOnCurve(Q)) throw PointNotOnCurveException(Q,EC);
    }
    catch (PointNotOnCurveException& e) {
        cout << "Caught an exception: " << e.what() << endl;
        return ECPoint(); // TODO : add something better than point at infinity here 
    }

    // cout << P.toString() << " + " << Q.toString() << " = "; 

    if (P.isPointAtInfinity) return Q; 
    if (Q.isPointAtInfinity) return P; 

    if (P.x == Q.x && (P.y != Q.y || P.y == ZZ(0)))
        return ECPoint(); // (P + (-P) = Identity point)

    if (P==Q){ // P = Q -->  P + Q = P + P 
        return doublePoint(P,EC);
    }

    ZZ_p::init(EC.p);

    ZZ_p P_x = to_ZZ_p(P.x);
    ZZ_p P_y = to_ZZ_p(P.y);

    ZZ_p Q_x = to_ZZ_p(Q.x);
    ZZ_p Q_y = to_ZZ_p(Q.y);

    ZZ_p lambdaNumerator = Q_y - P_y;
    ZZ_p lambdaDenominator = inv(Q_x - P_x);

    ZZ_p lambda = lambdaNumerator * lambdaDenominator;

    ZZ_p x_p = sqr(lambda) - P_x - Q_x ;
    ZZ x = rep(x_p);
    ZZ y = rep( lambda*(P_x - x_p) - P_y);

   return ECPoint(x,y);
}

// Scalar multiplication: kP
ECPoint scalarPointMultiplication(const ZZ& k_const, const ECPoint& P, const EllipticCurve& EC) {
    try {
        if (!EC.isOnCurve(P)) throw PointNotOnCurveException(P,EC);
    }
    catch (PointNotOnCurveException& e) {
        cout << "Caught an exception: " << e.what() << endl;
        return ECPoint(); // TODO : add something better than point at infinity here 
    }

    ECPoint result; 

    // done so we have modifiable lvalue for shifting k's bits and still be able to call function with ZZ(123)
    ZZ k = ZZ(k_const);  

    ECPoint base = P;
    turnOffCout();
    while (k > 0) {
        if (bit(k, 0)) { // check current LSB 
            result = addPoints(result, base, EC);
            cout << result.toString() << " | ";
        }
        base = doublePoint(base,EC);
        k >>= 1;
    }
    turnOnCout();
    // cout << endl <<k_const << "*" << P.toString() << " = " << result.toString() << ", ";
    return result;
}

ECPoint inversePoint(const ECPoint& P, const EllipticCurve& EC) {

    try {
        if (!EC.isOnCurve(P)) throw PointNotOnCurveException(P,EC);
    }
    catch (PointNotOnCurveException& e) {
        cout << "Caught an exception: " << e.what() << endl;
        return ECPoint(); // TODO : add something better than point at infinity here 
    }

    if (P == ECPoint()) return ECPoint();
    if (P.y == ZZ(0)) return P;

    ZZ_p y = to_ZZ_p(P.y);
    y = -y;
    return ECPoint(P.x,rep(y));   
}

ZZ bruteSearchSubgroupOrder(const ECPoint& P, const EllipticCurve& EC){
    if (!EC.isOnCurve(P)) {
        throw PointNotOnCurveException(P, EC);
        return (ZZ(0));
    }

    cout << "brute search order of generated subgroup from point " << P.toString() <<  endl;
    ECPoint temp = ECPoint();
    ZZ order = ZZ(0);
    // turnOffCout();
    do{
        temp = addPoints(temp,P,EC);
        order++;
    } while (!(temp == ECPoint()));
    // turnOnCout();
    return order;
}

ZZ bruteSearchPointMultiplicity(const ECPoint& P, const ECPoint& multipleOfP, const EllipticCurve& EC){

    //! WARNING :  user has to guarantee that P and multipleOfP is in the same subgroup, otherwise 
    //! the function will never end 

    if (!EC.isOnCurve(P)) throw PointNotOnCurveException(P, EC);
    if (!EC.isOnCurve(multipleOfP)) throw PointNotOnCurveException(P, EC);

    cout << "brute search point multiplicity of point " << P.toString() <<  endl;
    ECPoint temp = ECPoint();
    ZZ order = ZZ(0);
    // turnOffCout();
    do{
        temp = addPoints(temp,P,EC);
        order++;
        // cout << order << endl;
    } while (!(temp == multipleOfP));
    // turnOnCout();
    return order;
}

ZZ babyStepGiantStepECDLP(const ECPoint& P, const ECPoint& Q, const EllipticCurve& EC, const ZZ& subgroupOrder){
    // ZZ m = SqrRoot(subgroupOrder);
    ZZ m = SqrRoot(subgroupOrder) + ZZ(1);
    // cout << "number of baby steps to be computed: sqrt(" << subgroupOrder << ") :" << m << endl;
    unordered_map<ECPoint, ZZ, ECPointHasher> babySteps;

    // Compute baby steps: P, 2P, 3P, ..., jP
    ECPoint currentBabyStep = P;
    // cout << "Computing baby steps ..." << endl;
    // turnOffCout();

    // ZZ count = ZZ(0);
    for (ZZ j = ZZ(1); j <= m; ++j) {
        babySteps[currentBabyStep] = j;
        // cout << endl << " adding currentBabyStep : " << currentBabyStep.toString() << " with P : " << P.toString() << endl;
        currentBabyStep = addPoints(currentBabyStep, P, EC); 
        // count ++ ;
        // cout << "baby step counter : " << count << endl;
        // if (count == 10000){
        //     break;
        // }
    }
    // 1208925819614629174706176

    // Compute giant steps: Q, Q - mP, Q - 2mP, ..., Q - imP
    // * IDEA : precompute -mP beforehand
    // ECPoint currentGiantStep = Q;
    // cout <<  "Computing giant steps...";
    // turnOffCout();

    // ECPoint minusMP = inversePoint(scalarPointMultiplication(m,P,EC),EC);
    // for (ZZ i = ZZ(0); i <= m; ++i) {
    //     ECPoint curr = addPoints(Q,scalarPointMultiplication(i,minusMP,EC),EC);
    //     if (babySteps.find(curr) != babySteps.end()){
    //         ZZ j = babySteps[curr];
    //         ZZ k = i*m + j;
    //         return k;             
    //     }
    // }

    ECPoint minusMP = inversePoint(scalarPointMultiplication(m, P, EC), EC);
    ECPoint giantStepComponent = ECPoint();
    ECPoint currentGiantStepBase = Q; 
    for (ZZ i = ZZ(0); i <= m; ++i) {
        ECPoint curr = currentGiantStepBase; 

        if (babySteps.find(curr) != babySteps.end()){
            ZZ j = babySteps[curr];
            ZZ k = (i * m + j) % subgroupOrder;
            return k;              
        }

        currentGiantStepBase = addPoints(currentGiantStepBase, minusMP, EC); 
    }
    return ZZ(-1); 
}


struct PollardRhoState {
    ECPoint X;
    ZZ a, b;
};

ZZ partition(const ECPoint& P) {
    return P.x % ZZ(3); 
}

// ai bi P  == aj bj Q 
PollardRhoState f(const PollardRhoState& state, const ECPoint& P, const ECPoint& Q, const EllipticCurve& EC, const ZZ& subgroupOrder) {
    PollardRhoState nextState = state;
    ZZ category = partition(state.X);

    if (category == 0) {
        nextState.X = addPoints(state.X, P, EC);
        nextState.a = (state.a + 1) % subgroupOrder; 
    } else if (category == 1) {
        nextState.X = doublePoint(state.X,EC);
        nextState.a = (2 * state.a) % subgroupOrder; 
        nextState.b = (2 * state.b) % subgroupOrder; 
    } else {
        nextState.X = addPoints(state.X, Q, EC);
        nextState.b = (state.b + 1) % subgroupOrder; 
    }
    return nextState;
}

ZZ pollardRhoECDLP(const ECPoint& P, const ECPoint& Q, const EllipticCurve& EC, const ZZ& subgroupOrder) {

    PollardRhoState turtle = {P, ZZ(1), ZZ(0)}; 
    PollardRhoState rabbit = {P, ZZ(1), ZZ(0)}; 

    turtle = f(turtle, P, Q, EC, subgroupOrder);
    rabbit = f(f(rabbit, P, Q, EC, subgroupOrder), P, Q, EC, subgroupOrder);

    ZZ i = ZZ(0);
    while (!(turtle.X == rabbit.X)) {
        turtle = f(turtle, P, Q, EC, subgroupOrder);
        rabbit = f(f(rabbit, P, Q, EC, subgroupOrder), P, Q, EC, subgroupOrder);
        i++;
    }
         
    // cout << endl << "Collision found at iteration: " << i << endl;
    // cout << "Turtle X: " << turtle.X.toString() << ", a=" << turtle.a << ", b=" << turtle.b << endl;
    // cout << "Rabbit X: " << rabbit.X.toString() << ", a=" << rabbit.a << ", b=" << rabbit.b << endl;

    ZZ numerator = (turtle.a - rabbit.a) % subgroupOrder;
    ZZ denominator = (rabbit.b - turtle.b) % subgroupOrder;

    if (denominator == 0) {
        // cout << "Denominator is zero, no inverse. Try again or different f function." << endl;
        return ZZ(-1);
    }

    ZZ_p::init(subgroupOrder); 


    ZZ_p p_inv = inv(to_ZZ_p(denominator));

    if (p_inv == 0) {
        cout << "Denominator has no inverse modulo subgroupOrder." << endl;
        return ZZ(-1);
    }
    // cout << "Inverse of denominator mod subgroup order:  " << rep(p_inv) << endl;

    ZZ k = rep((to_ZZ_p(numerator) * p_inv));

    cout << "Pollard Rho found potential k: " << k << endl;
    return k;
        
}

void solveECDLP_BothMethods(const EllipticCurve& curve, const ECPoint& P, const ZZ& k,const ZZ& order) {


    cout << "_____SOLVING ECDLP WITH ALL METHODS_____" << endl;
    cout << curve.toString() << endl;

    // try {
    //     if (!curve.isOnCurve(P)) throw PointNotOnCurveException(P,curve);
    // }
    // catch (PointNotOnCurveException& e) {
    //     cout << "Caught an exception: " << e.what() << endl;
    //     return; 
    // }
    
    // curve.print();
    cout << "P : " << P.toString() << endl;
    cout << "k : " << k << endl;

    ECPoint Q = scalarPointMultiplication(k, P, curve);

    // cout << " k : " << k << "  P: " <<  P.toString() << " curve : " << curve.toString() << endl;
    cout << endl << endl <<  "Q : " << Q.toString() << endl;

    cout << "k*P=" << Q.toString() << endl;
    //* brute search subgroup order (we would ideally know it upfront) 
    ZZ subgroupOrder = order;
    if (order == ZZ(0)){
        cout <<"we are brute searching subgroup order because it was not specified" << endl;
        subgroupOrder = bruteSearchSubgroupOrder(P, curve);
        cout << "done brute searching subgroup order" << endl;
    }

    cout << "k % order = " << k << " % " << subgroupOrder << " = " << k%subgroupOrder << endl;

    //* brute search multiplicity of P (should be same as k) 
    // ZZ multiplicity = bruteSearchPointMultiplicity(P,Q,curve);
    // cout << "multiplicity of P :" << multiplicity << endl;
    // cout << "brute search found correct k : " << ((multiplicity==k%subgroupOrder)? "Yes" : "No") << endl;

    cout << endl << "--- Solving ECDLP using Baby Step Giant Step ---" << endl;
    ZZ BSGS_k = babyStepGiantStepECDLP(P, Q, curve, subgroupOrder);
    cout << endl << "Baby Step Giant Step result (k): " << BSGS_k << endl;
    ECPoint BSGSMultiple = scalarPointMultiplication(BSGS_k, P, curve);
    cout << "Verification (BSGS result * P == Q): " << ((BSGSMultiple == Q) ? "Yes" : "No") << endl;

    cout << endl << "--- Solving ECDLP using Pollard Rho ---" << endl;
    ZZ PollardRho_k = pollardRhoECDLP(P, Q, curve, subgroupOrder);
    cout << endl << "Pollard Rho result (k): " << PollardRho_k << endl;
    ECPoint pollardRhoMultiple = scalarPointMultiplication(PollardRho_k, P, curve);
    cout << "Verification (Pollard Rho result * P == Q): " << ((pollardRhoMultiple == Q) ? "Yes" : "No") << endl;
}



// void printBitsMSBtoLSB(const ZZ& num) {
//     long numBits = NumBits(num); 

//     cout << "Bits from MSB to LSB: ";
//     for (long i = numBits - 1; i >= 0; i--) { 
//         cout << bit(num, i); 
//     }
//     cout << endl;
// }

void printZZ(const ZZ& x) {
    cout << x << endl;
}


// Function to extract bits from a ZZ number and store them in a vector<bool>
vector<bool> getBitVector(const ZZ& num, long bitSize) {
    vector<bool> bitVec(bitSize, false); // Initialize with false (0s)
    
    long numBits = NumBits(num); // Get actual bit length of the number
    for (long i = 0; i < numBits; i++) {
        bitVec[bitSize - 1 - i] = bit(num, i); // Store bit at the correct index
    }
    
    return bitVec;
}

// Function to align two numbers into bit vectors with equal length
pair<vector<bool>, vector<bool>> alignNumbers(const ZZ& a, const ZZ& b) {
    long maxBits = max(NumBits(a), NumBits(b)); // Determine required bit length
    
    vector<bool> aBits = getBitVector(a, maxBits); // Get aligned bits
    vector<bool> bBits = getBitVector(b, maxBits);

    return {aBits, bBits}; // Return aligned vectors
}

// computing s*P + t*Q
ECPoint shamirsTrick(const ECPoint& P, const ECPoint& Q, const ZZ& s, const ZZ& t, const EllipticCurve& EC) {
    auto [sBits, tBits] = alignNumbers(s, t);
    long bitLength = sBits.size();
    
    ECPoint result = ECPoint();
    ECPoint sumOfPandQ = addPoints(P,Q,EC);

    for (long i = 0; i < bitLength; i++) {
        result = doublePoint(result,EC);

        if (sBits[i] && tBits[i]) {
            result = addPoints(result,sumOfPandQ,EC);
        } else if (sBits[i]) {
            result = addPoints(result,P,EC);
        } else if (tBits[i]) {
            result = addPoints(result,Q,EC);
        }
    }
    
    return result;
}

// Function to print bit vectors
void printBitVector(const vector<bool>& bits) {
    for (bool bit : bits) {
        cout << bit;
    }
    cout << endl;
}

struct Curve {
    ZZ a;
    ZZ b;
};

struct Generator {
    ZZ x;
    ZZ y;
};

struct ECCurveData {
    int bit_size;
    ZZ prime;
    Curve curve;
    ZZ curve_order;
    Generator generator;
    ZZ generator_order;
    int log2_order;
};


#include <iostream>
#include <fstream>
#include <vector>
#include <nlohmann/json.hpp>
#include <NTL/ZZ.h>

using json = nlohmann::json;

struct CurveParams {
    long bit_size;
    ZZ prime;
    ZZ curve_a;
    ZZ curve_b;
    ZZ curve_order;
    ZZ generator_x;
    ZZ generator_y;
    ZZ generator_order;
    long log2_order;
};

int main() {
        
    // ZZ a = ZZ(40);
    // ZZ b = ZZ(84);
    // ZZ p = ZZ(97);
    
    // EllipticCurve curve(a, b, p);

    // curve.print();

    // ECPoint infinityPoint = ECPoint();

    // //! TEST POINTS
    // // Define points P and Q ( P = Q^(-1) ... if connected P and Q form a vertical line)

    // ECPoint P = ECPoint(ZZ(8), ZZ(25)); 
    // ECPoint Q = ECPoint(ZZ(8), ZZ(72)); 
    
    // // Edge point is any point inverse to itself ( y coord = 0)
    // ECPoint EdgePoint = ECPoint(ZZ(21), ZZ(0));
        
    // // Define points R and S ( R != S^(-1) ... if connected P and Q do NOT form a vertical line)
    // ECPoint R = ECPoint(ZZ(9), ZZ(94)); 
    // ECPoint S = ECPoint(ZZ(11), ZZ(20));

    // cout << "Point P: " << P.toString();
    // cout << "Point Q: " << Q.toString();

    // cout << "Is P on curve? " << (curve.isOnCurve(P) ? "Yes" : "No") << endl;
    // cout << "Is Q on curve? " << (curve.isOnCurve(Q) ? "Yes" : "No") << endl;

    // //! BASIC EC ARITHMETICS AND FUNCTIONS
    // // 0 + P 
    // cout << infinityPoint.toString() << " + " << P.toString() << " = "<< addPoints(infinityPoint,P,curve).toString() << endl;
    // // P + 0
    // cout << P.toString() << " + " << infinityPoint.toString() << " = "<< addPoints(P,infinityPoint,curve).toString() << endl;
    // // 0 + 0
    // cout << infinityPoint.toString() << " + " << infinityPoint.toString() << " = "<< addPoints(infinityPoint,infinityPoint,curve).toString() << endl;
    
    // // Q + P ... (P = inverse(Q))
    // cout << Q.toString() << " + " << P.toString() << " = "<< addPoints(Q,P,curve).toString() << endl;
    // // P + Q ...(P = inverse(Q))
    // cout << P.toString() << " + " << Q.toString() << " = "<< addPoints(P,Q,curve).toString() << endl;
    
    // // P + Q ... (P != inverse(Q)) 
    // cout << R.toString() << " + " << S.toString() << " = "<< addPoints(R,S,curve).toString() << endl;
    // // Q + P ... (P != inverse(Q)) 
    // cout << S.toString() << " + " << R.toString() << " = "<< addPoints(S,R,curve).toString() << endl;

    // // P + P(^-1) ... (P = inverse(P))  
    // cout << addPoints(EdgePoint,EdgePoint,curve).toString() << endl;
    
    // ZZ order = bruteSearchSubgroupOrder(R,curve);

    // cout << "order of subgroup generated by "<< R.toString() << " is " << order << endl;

    // // (additive) inverse point calculation
    // cout << inversePoint(S,curve).toString() << endl;
    // cout << inversePoint(EdgePoint,curve).toString() << endl;
    // cout << inversePoint(infinityPoint,curve).toString() << endl;
    
    // // override of == operator tests
    // cout << ECPoint().toString() << " == " << infinityPoint.toString() << " --> " <<(ECPoint() == infinityPoint ? "Yes" : "No") << endl;
    // cout << infinityPoint.toString() << " == " << infinityPoint.toString() << " --> " <<(infinityPoint == infinityPoint ? "Yes" : "No") << endl;
    // cout << R.toString() << " == " << infinityPoint.toString() << " --> " <<(R == infinityPoint ? "Yes" : "No") << endl;
    // cout << ECPoint(ZZ(0),ZZ(0)).toString() << " == " << infinityPoint.toString() << " --> " << (ECPoint(ZZ(0),ZZ(0)) == infinityPoint ? "Yes" : "No") << endl;
    
    // //! HASHER FUNCTIONALITY
    // ECPointHasher hasher;

    // ECPoint p3= ECPoint(ZZ(1), ZZ(2));
    // ECPoint p4= ECPoint(ZZ(1), ZZ(2));
    // ECPoint p5= ECPoint(ZZ(2), ZZ(2));
    
    // cout << "hash("<< p3.toString()<<") == hash("<<p4.toString()<<") ? "<<((hasher(p3) == hasher(p4))? "Yes" : "No") << endl;  
    // cout << "hash("<< p5.toString()<<") == hash("<<p4.toString()<<") ? "<<((hasher(p5) == hasher(p4))? "Yes" : "No") << endl;  
    
    // cout << ((ECPointHasher{}(ECPoint(ZZ(1),ZZ(2))) == ECPointHasher{}(ECPoint(ZZ(1),ZZ(2))))? "True" : "False") << endl; 
    
    // unordered_map<ECPoint, ZZ, ECPointHasher> pointMap;

    // ECPoint X = ECPoint(ZZ(10), ZZ(20));  
    // ECPoint Z= ECPoint(ZZ(15), ZZ(25));  

    // pointMap[X] = ZZ(5);  
    // pointMap[Z] = ZZ(10); 

    // cout << "Hash of " << X.toString() << " : "<<  ECPointHasher{}(X) << endl;
    // cout << "Hash of " << Z.toString() << " : "<< ECPointHasher{}(Z) << endl;

    // cout << "Value associated with point X: " << pointMap[X] << endl;
    // cout << "Value associated with point Z: " << pointMap[Z] << endl;

    // // auto start = std::chrono::high_resolution_clock::now();
    // // <CODE FOR TIME TESTING HERE>
    // // auto end = std::chrono::high_resolution_clock::now();
    // // std::chrono::duration<double> duration = end - start;    
    // // std::cout << "Execution time: " << duration.count() << " seconds" << endl;
    
    // //! CURVE WITH NON-PRIME ORDER SUBGROUP
    // EllipticCurve curve_nonPrimeOrderSubgroup(ZZ(1), ZZ(1), ZZ(997));    
    // ECPoint P_nonPrimeOrderSubgroup(ZZ(278), ZZ(599)); // order of subgroup generated by P is NOT prime ...
    // ZZ k_nonPrimeOrderSubgroup = to_ZZ("3012385230498523049528334985793487132487149879238472983742983410000239874501");
    
    // solveECDLP_BothMethods(curve_nonPrimeOrderSubgroup,P_nonPrimeOrderSubgroup,k_nonPrimeOrderSubgroup,ZZ(995));
    
    // //! CURVE WITH PRIME ORDER SUBGROUP
    // EllipticCurve curve_primeOrderSubgroup(ZZ(1), ZZ(2), ZZ(827));
    // ECPoint P_primeOrderSubgroup(ZZ(774), ZZ(437));
    // ZZ k_primeOrderSubgroup = ZZ(155523);
    // solveECDLP_BothMethods(curve_primeOrderSubgroup,P_primeOrderSubgroup,k_primeOrderSubgroup,ZZ(107));
    
    
    // // for (ZZ i = ZZ(0); i < to_ZZ("20000000000000000000000000000"); i+= to_ZZ("100000000000000000000"))
    // // {
    // //     cout << " ahoj " << i <<  endl;

    // //     if(i==to_ZZ("10000000000000000000000000")){
    // //         cout << "it has happened" << i << endl;
    // //     }
    // // }
    
    // for (int i = 0; i < 3; i++)
    // {
    //     cout << " TEST NORMAL INT " << i <<  endl;
    //     continue;
    // }
    
    // for (ZZ i = ZZ(0); i < ZZ(3); i++)
    // {
    //     cout << " TEST ZZ NUM " << i <<  endl;
    //     continue;
    // }

    
    // EllipticCurve curve_primeOrderSubgroup2(to_ZZ("1461501637330902918203684832716283019653785059324"), to_ZZ("163235791306168110546604919403271579530548345413"), to_ZZ("1461501637330902918203684832716283019653785059327"));
    // EllipticCurve curve_primeOrderSubgroup2(ZZ(2), ZZ(3), ZZ(100003));
    // ECPoint P_primeOrderSubgroup2(to_ZZ("11549"), to_ZZ("23288"));
    // // ZZ k_primeOrderSubgroup2 = ZZ(15523);
    // solveECDLP_BothMethods(curve_primeOrderSubgroup2,P_primeOrderSubgroup2,ZZ(15523),ZZ(100294));
    
    // //! CURVE WITH PRIME ORDER SUBGROUP
    // EllipticCurve curve_nonPrimeOrderSubgroup2 = EllipticCurve(to_ZZ("1461501637330902918203684832716283019653785059324"), to_ZZ("163235791306168110546604919403271579530548345413"), to_ZZ("1461501637330902918203684832716283019653785059327"));    
    // ECPoint P_nonPrimeOrderSubgroup2(to_ZZ("425826231723888350446541592701409065913635568770"), to_ZZ("203520114162904107873991457957346892027982641970")); // order of subgroup generated by P is NOT prime ...
    // // ZZ k_nonPrimeOrderSubgroup2 = to_ZZ("529843787465298374652983746592873465982374619823");
    // ZZ k_nonPrimeOrderSubgroup2 = to_ZZ(1);

    // ZZ order2 = to_ZZ("1461501637330902918203687197606826779884643492439");

    std::ifstream file("testcurves.json");
    if (!file) {
        std::cerr << "Failed to open JSON file.\n";
        return 1;
    }

    json j;
    try {
        file >> j;
    } catch (const json::parse_error& e) {
        std::cerr << "JSON parse error: " << e.what() << '\n';
        return 1;
    }

    if (!j.is_array()) {
        std::cerr << " Expected a JSON array at the top level.\n";
        return 1;
    }

    vector<CurveParams> curves;

    try {
        for (const auto& entry : j) {
            CurveParams cp;

            cp.bit_size = entry["bit_size"].get<long>();
            cp.prime = ZZ(entry["prime"].get<long>());
            cp.curve_a = ZZ(entry["curve"]["a"].get<long>());
            cp.curve_b = ZZ(entry["curve"]["b"].get<long>());
            cp.curve_order = ZZ(entry["curve_order"].get<long>());
            cp.generator_x = ZZ(entry["generator"]["x"].get<long>());
            cp.generator_y = ZZ(entry["generator"]["y"].get<long>());
            cp.generator_order = ZZ(entry["generator_order"].get<long>());
            cp.log2_order = entry["log2_order"].get<long>();

            curves.push_back(cp);
        }
    } catch (const json::type_error& e) {
        std::cerr << " Type error while reading JSON: " << e.what() << '\n';
        return 1;
    }

    clock_t start;
    clock_t end;
    clock_t startbrute;
    clock_t endbrute;


    vector<double> rho_results;
    vector<double> bsgs_results;
    vector<double> bf_results;
    vector<double> bruteforce_results;

    std::vector<std::string> values;

    std::stringstream rho_ss;
    std::stringstream bsgs_ss;
    std::stringstream bruteforce_ss;

    string rho_commaSeparatedString;
    string bsgs_commaSeparatedString;
    string bruteforce_commaSeparatedString;

    std::ofstream outputFile("new_output.csv", std::ios::app);

    for (size_t i = 0; i < curves.size(); ++i) {
        const auto& curve = curves[i];
        std::cout << "Curve #" << i + 1 << ":\n"
                  << "  Bit Size:        " << curve.bit_size << "\n"
                  << "  Prime:           " << curve.prime << "\n"
                  << "  Curve a:         " << curve.curve_a << "\n"
                  << "  Curve b:         " << curve.curve_b << "\n"
                  << "  Curve Order:     " << curve.curve_order << "\n"
                  << "  Generator X:     " << curve.generator_x << "\n"
                  << "  Generator Y:     " << curve.generator_y << "\n"
                  << "  Generator Order: " << curve.generator_order << "\n"
                  << "  Log2 Order:      " << curve.log2_order << "\n\n";

        ECPoint Point = ECPoint(curve.generator_x,curve.generator_y);
        EllipticCurve EC = EllipticCurve(curve.curve_a,curve.curve_b,curve.prime);
        cout << "key Bit Size: " << curve.bit_size << endl;

        ZZ privateKey = NTL::RandomBnd(curve.curve_order) + 1;
        cout << "curve order : " << curve.curve_order << endl;
        ECPoint multipliedPoint = scalarPointMultiplication(privateKey,Point,EC);
    
        start = clock();
        ZZ result = pollardRhoECDLP(Point,multipliedPoint,EC,curve.curve_order);
        // ZZ result = bruteSearchPointMultiplicity(Point,multipliedPoint,EC);
        end = clock();
        cout << "privateKey: " << privateKey << endl;
        clock_t elapsed = end - start;        
        bruteforce_results.push_back(elapsed);
        // rho_results.push_back(elapsed);
        // printf("Time measured by Pollard's Rho method: %ld seconds.\n", elapsed);
        printf("Time measured by brute force method: %ld seconds.\n", elapsed);
        cout << "Brute force result " << result << endl; 
        cout << ((privateKey % curve.prime == result % curve.prime)? "GOOD" : "NOT FOUND" ) << endl << endl; 
    
        // start = clock();
        // result = babyStepGiantStepECDLP(Point,multipliedPoint,EC,curve.curve_order);
        // end = clock();
        // // elapsed = double(end - start)/CLOCKS_PER_SEC;
        // elapsed = end-start;
        // bsgs_results.push_back(elapsed);        
        // printf("Time measured by BSGS method: %ld .\n", elapsed);
        // cout << "BSGS method result " << result << endl; 
        // cout << ((privateKey % curve.prime == result % curve.prime)? "GOOD" : "NOT FOUND" ) << endl << endl; 

        // rho_ss << int(rho_results[i]);
        // bsgs_ss << int(bsgs_results[i]);
        bruteforce_ss << int(bruteforce_results[i]);

        if(i%100 != 99){
            // rho_ss << ",";
            // bsgs_ss << ",";
            bruteforce_ss << ",";
        }
        else{
            // rho_commaSeparatedString = rho_ss.str();
            // bsgs_commaSeparatedString = bsgs_ss.str();
    
            // std::ofstream outputFile("new_output.csv", std::ios::app);
    
            // if (outputFile.is_open()) {
                // outputFile << rho_ss.str() << "\n";
                // outputFile << bsgs_ss.str() << "\n";
            outputFile << bruteforce_ss.str() << "\n";
                
                // std::cout << "String appended to 'output.csv'" << std::endl;
            // } else {
            //     std::cerr << "Unable to open file for writing." << std::endl;
            // }

            rho_ss = stringstream();
            bsgs_ss = stringstream();
            bruteforce_ss = stringstream();
        }

    }

    outputFile.close();

    // std::stringstream rho_ss;
    // std::stringstream bsgs_ss;
    

    // for (size_t i = 0; i < rho_results.size(); ++i) {
    //     rho_ss << int(rho_results[i]);
    //     bsgs_ss << int(bsgs_results[i]);
    //     if (i < rho_results.size() - 1) {
    //         rho_ss << ",";
    //         bsgs_ss << ",";
    //     }
    // }

    // string rho_commaSeparatedString = rho_ss.str();
    // string bsgs_commaSeparatedString = bsgs_ss.str();

    // std::cout << "Comma-separated string for rho: " << rho_commaSeparatedString << std::endl;
    // std::cout << "Comma-separated string for bsgs: " << bsgs_commaSeparatedString << std::endl;


    // if (outputFile.is_open()) {
    //     outputFile << rho_commaSeparatedString << "\n";
    //     outputFile << bsgs_commaSeparatedString << "\n";
    //     outputFile.close();
    //     std::cout << "String appended to 'output.csv'" << std::endl;
    // } else {
    //     std::cerr << "Unable to open file for writing." << std::endl;
    // }


    double rho_avg=0.0;
    double bsgs_avg=0.0;

    for(double entry : rho_results){
        rho_avg += entry;
    }

    rho_avg = rho_avg / 100;

    for(double entry : bsgs_results){
        bsgs_avg += entry;
        
    }

    bsgs_avg = bsgs_avg / 100;


    cout << "rho average : " << rho_avg << endl;
    cout << "bsgs average : " << bsgs_avg << endl;
    return 0;

    

    // pani sabikova zabazpeci podpisy riaditela atd 

    // for (const auto& curve : curves) {
        // ECPoint Point = ECPoint(curve.generator.x,curve.generator.y);
        // EllipticCurve EC = EllipticCurve(curve.curve.a,curve.curve.b,curve.prime);


        // cout << "Prime: " << curve.prime << endl;
        // cout << "Curve a: " << curve.curve.a << ", b: " << curve.curve.b << endl;
        // cout << "Curve Order: " << curve.curve_order << endl;
        // cout << "Generator x: " << curve.generator.x << ", y: " << curve.generator.y << endl;
        // cout << "Generator Order: " << curve.generator_order << endl;
        // cout << "Log2 Order: " << curve.log2_order << endl;
        // cout << "---------------------------------" << endl;

        // ZZ privateKey = curve.curve_order - ZZ(200);


        

        // start = clock();

        // ZZ result = bruteSearchPointMultiplicity(Point,multipliedPoint,EC);

        // end = clock();
        
        // double elapsedbrute = double(end - start)/CLOCKS_PER_SEC;        
        // printf("Time measured by brute force method: %.3f seconds.\n", elapsedbrute);
        // cout << "private key " << privateKey << endl;

        // cout << "brute force result " << result << endl; 
        // cout << ((privateKey % curve.prime == result % curve.prime)? "GOOD" : "NOT FOUND" ) << endl; 

        
        // bruteSearchPointMultiplicity()

    // }


    // double sum = 0;
    // double add = 1;

    // // Start measuring time
    // start = clock();
    
    // // int iterations = 1000*1000*1000;
    // // for (int i=0; i<iterations; i++) {
    // //     sum += add;
    // //     add /= 2.0;
    // // }

    // // Stop measuring time and calculate the elapsed time
    // end = clock();
    // double elapsed = double(end - start)/CLOCKS_PER_SEC;
    
    // printf("Result: %.20f\n", sum);
    
    // printf("Time measured: %.3f seconds.\n", elapsed);
    
    // return 0;
}


// stack exchange parameters 
//Let E:y^2=x^3+8x+7 over F_73. There are 82 points on this curve. So the only possible subgroup sizes are 2 and 41. Sure enough, if P=(32,53), you have [41]P=O, and 41|82.