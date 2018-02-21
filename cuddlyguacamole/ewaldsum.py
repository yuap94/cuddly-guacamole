


def k_energy (i, r, q, k_c):
    
    #prefactor 
    V=np.prod(box)
    pre_fac = 1/(2*V*epsilon_0)  
    
    for i in range (n):  
        for j in range (m): 
            
            
            for k_i in range (-k_c,k_c+1)
                for k_j in range (-k_c,k_c+1)
                    for k_k in range (-k_c,k_c+1) 
                        if np.linalg.norm(k_i,k_j,k_k) <= k_c 
                        #Reciprocal vector 
                        k = 2.0*np.pi*np.array([k_i / (box[0]), k_j / (box[1]), k_k / (box[2])] )                               
                        k_sq = k*k  
    
                        exp_term += np.exp(-(sigma**2)*k_sq/2)/k_sq) 
     
                k_space_energy= pre_fac*(s_k**2)*exp_term
    
    return k_space_energy


def k_cut_off (self) :
    #!!! accuracy set to 1e-6 gives p, make p global 
    k_c = (2*p)/r_cut 
    
    return k_c 

def sigma (self) : 
 
    sigma = r_c / np.sqrt (2*p) 
    
    return sigma 
    
                     
def str_fac (self) :
    
    for i in range ():
        s_k=0 
        
        for j in range (p_num) :
            s_k += charge(j)*np.cos(np.dot(k,position(j)))
            
    return s_k 


