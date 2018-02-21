


def Ewald_long_energy (i, r, q, k_c):
    
    #prefactor 
    V=np.prod(box)
    pre_fac = 1/(2*V*epsilon_0)  
    
    long_energy = 0 
    str_factor = 0 
    for i in range (len(k_vector)):  
        k=k_vector[i]
        #length of k
        k_length = np.linalg.norm (k) 
        k_length2 = k_length**2 
        
        for j in range (0,particle_number):
             q = charge[i]
             r = coord [i] 
             str_fac += q*np.cos(np.dot(k,r))
             exp_term += np.exp(-(sigma**2)*k_length2/2)/k_length2 
     
        long_energy += np.linalg.norm(str_fac)**2*exp_term 
    long_energy *= pre_fac 
    
    return long_energy

def k_vectors(): 
k_vector=[]
   for k_i in range (-k_c,k_c+1):
       for k_j in range (-k_c,k_c+1):
           for k_k in range (-k_c,k_c+1):
                if np.linalg.norm(k_i,k_j,k_k) <= k_c 
                #Reciprocal  lattice vector 
                k = 2.0*np.pi*np.array([k_i / (box[0]), k_j / (box[1]), k_k / (box[2])] )   
                k_vector.append(k) 
                        
                        
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


