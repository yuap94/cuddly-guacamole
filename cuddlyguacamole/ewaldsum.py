


def k_energy (i, r, q, cutoff_kspace):
    
    
    for k_i in range (-k_c,k_c+1)
        for k_j in range (-k_c,k_c+1)
            for k_k in range (-k_c,k_c+1) 
                #Reciprocal vector 
                k = 2.0*np.pi*np.array([k_i / (box[0]), k_j / (box[1]), k_k / (box[2])] )                               
                k2 = k*k  
    
    #prefactor, V value ???
    pre_fac = 1/(2*V*epsilon_0)  

    #Structure factor 
    
    s_k = sum( q*np.exp(k*r))          
    s_k2= s_k**2              

    exp_term = np.sum(exp(-sigma2*k2/2)/k2) 
     
         k_space_energy= pre_fac*s_k2*exp_term
        
    return k_space_energy


def k_cut_off (self) :
    #accuracy set to 1e-6 gives p 
    p = ????
    k_c = 2*p/r_cut 
    #Sigma
    sigma = r_c / np.sqrt (2*p) 
    sigma2 = sigma*sigma 
    

                      
 
def str_fac (self) :
    
    for i in range ():
        s_k=0 
        
        for j in range (p_num) :
            s_k += charge(j)*np.cos(np.dot(k,position(j)))
            
    return s_k 


