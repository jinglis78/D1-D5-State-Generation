from collections import defaultdict

# chi(j) = sum^(m=j)_(m=-j) y^(2m)
# output = dicotionary with key = y coeff, value = number of terms at that order
    
def chi(s, N=0, y_power=1):
    terms = {}
    # multiply everything by 2 to avoid fractions
    for m2 in range(int(-2*s), int(2*s) + 1, 2):  
        m = m2 / 2.0
        exponent = ( 2 * m + N ) * y_power 
        coefficient = 1
        if exponent in terms:
            terms[exponent] += coefficient
        else:
            terms[exponent] = coefficient
    return terms

def chiTilde(s, N=0, y_power=1, y_val=1):
    """
    Compute χ_{s;N}(y) evaluated at y=y_val, treating it as a polynomial.
    
    Parameters:
        s       -- spin (can be half-integer, e.g. 0.5, 1, 1.5)
        N       -- shift exponent
        y_power -- multiplier for exponents
        y_val   -- numerical value of y (default = 1)
    
    Returns:
        value at y = y_val
    """
    value = 0
    for m2 in range(int(-2*s), int(2*s) + 1, 2):  
        m = m2 / 2.0
        exponent = (2*m + N) * y_power
        value += y_val**exponent
    return value

    
print (chi(1))
# phi_short(j) = q^j ( sum^(2j)_(m=0)  q^(m/2)  chi(j-n/2) ) / (1-q) 
# output = dictionary of dictionaries, where key1 = q_order, key2 = y_order, and value = coeff at that order

def phi_short(j, q_order, y_power = 1):
    # dictionary to store the terms
    terms = {}
    # q^j / (1-q) = q^j + q^(j+1) + q^(j+2) + ...
    """for i in range(0, q_order+1):
        if (i+j) < q_order:
            terms[i+j] = 1
    
    print(terms)"""
    
    sum_terms = {}
    max = int(2*j + 1)
    for m in range(0, max):  # inclusive of 2j
        chi_m = chi(j-m/2)
        sum_terms[m/2] = chi_m

    #print(sum_terms)
    
    for i in range(0, q_order+1):
        if (i+j) < q_order:
            #print (i+j)
            for k in sum_terms:
                # sum_terms[k] = y-expansion at order q^k
                
                for l in sum_terms[k]:
                    # l = y_coeff
                    # sum_terms[k][l] = coeff at order q^k qnd y^l
                    if ((k+i+j)*y_power) in terms:
                        if (l*y_power) in terms[(k+i+j)*y_power]:
                            terms[(k+i+j)*y_power][l*y_power] += sum_terms[k][l]
                        else:
                            terms[(k+i+j)*y_power][l*y_power] = sum_terms[k][l]
                    else:
                        terms[(k+i+j)*y_power] = {}
                        terms[(k+i+j)*y_power][l*y_power] = sum_terms[k][l]
                    
    

        
    return terms

def add_nested_dicts(d1, d2):
    result = defaultdict(lambda: defaultdict(int))
    
    # add everything from d1
    for k1, subdict in d1.items():
        for k2, val in subdict.items():
            result[k1][k2] += val
    
    # add everything from d2
    for k1, subdict in d2.items():
        for k2, val in subdict.items():
            result[k1][k2] += val
    
    # convert back to normal dict
    return {k1: dict(subdict) for k1, subdict in result.items()}

def multiply_polynomials(d1, d2):
    result = defaultdict(int)
    for exp1, coeff1 in d1.items():
        for exp2, coeff2 in d2.items():
            new_exp = exp1 + exp2
            new_coeff = coeff1 * coeff2
            result[new_exp] += new_coeff
    return dict(result)

def multiply_nested_polynomials(p1, p2):
    result = defaultdict(lambda: defaultdict(int))
    
    for q1, inner1 in p1.items():
        for q2, inner2 in p2.items():
            q_exp = q1 + q2
            # multiply inner dicts (the y-parts)
            for y1, coeff1 in inner1.items():
                for y2, coeff2 in inner2.items():
                    y_exp = y1 + y2
                    result[q_exp][y_exp] += coeff1 * coeff2
    
    # convert defaultdicts back to dicts
    return {q: dict(inner) for q, inner in result.items()}

def X(k, y_power = 1):
    phi1 = phi_short((k-1)/2, 4, y_power)
    phi2 = phi_short((k)/2, 4, y_power)
    phi3 = phi_short((k+1)/2, 4, y_power)
    
    chi1 = chiTilde(1/2, k, y_power)
    chi2 = chiTilde(0, k, y_power)
    
    phi1_mult = {}
    phi2_mult = {}
    phi3_mult = {}
    
    for itr in phi3:
            # phi3[k] = y-expansion at order q^k
            # multiply with chi2
            
            product = {itr1: itr2 * chi1 for itr1,itr2 in phi3[itr].items()}
            
            if (itr) in phi3_mult:
                phi3_mult[itr] = add_nested_dicts(phi3_mult[itr], product)
            else:
                phi3_mult[itr] = product
                
    for itr in phi2:
            # phi3[k] = y-expansion at order q^k
            # multiply with chi2
            
            product = {itr1: itr2 * chi2 for itr1,itr2 in phi2[itr].items()}
            
            if (itr) in phi2_mult:
                phi2_mult[itr] = add_nested_dicts(phi2_mult[itr], product)
            else:
                phi2_mult[itr] = product
    
    for itr in phi1:
            # phi3[k] = y-expansion at order q^k
            # multiply with chi2
            
            product = {itr1: itr2 * chi1 for itr1,itr2 in phi1[itr].items()}
            
            if (itr) in phi1_mult:
                phi1_mult[itr] = add_nested_dicts(phi1_mult[itr], product)
            else:
                phi1_mult[itr] = product
    
    sum_dict = add_nested_dicts(phi1_mult, phi2_mult)
    
    sum_dict = add_nested_dicts(sum_dict, phi3_mult)
    
    print(f"This is X({k})")
    return sum_dict

def pretty_print_expansion(X, q_power_max = 4):
    for q_power, y_terms in sorted(X.items(), key=lambda x: x[0]):
        # Build the y-polynomial as a list of strings
        if q_power_max > q_power:
            terms = []
            for y_power, coeff in sorted(y_terms.items(), key=lambda x: x[0]):
                if coeff == 0:
                    continue
                # format coefficient
                if coeff == 1:
                    coeff_str = ""
                else:
                    coeff_str = f"{coeff}*"
                # format y term
                if y_power == 0:
                    y_str = "1"
                elif y_power == 1:
                    y_str = "y"
                else:
                    y_str = f"y^{int(y_power) if y_power.is_integer() else y_power}"
                terms.append(f"{coeff_str}{y_str}")
            
            # Join terms with " + "
            y_poly = " + ".join(terms) if terms else "0"
            print(f"q^{q_power} * ({y_poly})")

print ("phi_short")
print(phi_short(3/2, 3))

X1 = X(2)


def X_squared(k):
    Xk = X(k)
    squared = multiply_nested_polynomials(Xk, Xk)
    return squared

def X_qsquared(k):
    Xk = X(k, y_power=2)
    return Xk

Xsquare = X_squared(2)
Xqsquare = X_qsquared(2)


print("This is X2:")
pretty_print_expansion(X1)

print("This is Xsquare:")
pretty_print_expansion(Xsquare)

print("This is Xqsquare:")
pretty_print_expansion(Xqsquare)

factor = 0.5
for k in Xqsquare:
    for k1 in  Xqsquare[k]:
        Xqsquare[k][k1] *= factor
        
for k in Xsquare:
    for k1 in  Xsquare[k]:
        Xsquare[k][k1]  *= factor


Xsum = add_nested_dicts(add_nested_dicts(Xqsquare, Xsquare),X1)

print(Xsum)

pretty_print_expansion(Xsum)
