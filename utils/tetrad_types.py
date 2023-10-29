TETRAD_MAP_FOUR = {
    ## Non- Crossover ##
    ('blue','blue','orange','orange'): 'A',
    ##########################
    ### Single Crossovers ####
    ##########################
    ('blue','orange','purple','yellow'): 'B',  
    ('blue','empty','orange', 'white'): 'C',
    ###########################
    #### Double Crossovers ####
    ###########################
    ('blue','green', 'orange', 'red'): 'D', 
    ('empty','green', 'orange', 'purple'): 'E',
    ('blue','red', 'white', 'yellow'): 'F',
    ('empty','purple', 'white', 'yellow'): 'G',
    ('purple','purple', 'yellow', 'yellow'): 'H',    
    ###########################
    #### Triple Crossovers ####
    ###########################
    ('green','purple','red','yellow') : 'I',
    ('empty','green','red','white') : 'J',
    ###########################
    #### Quadra Crossovers ####
    ###########################
    ('green','green', 'red', 'red') : 'K',
    ###########################
    ######### MI and MII NDJ ##########
    ###########################
    ('empty','empty', 'white', 'white'): 'L',
    ('blue','blue', 'empty', 'orange') : 'M', 
    ###############################
    ######  Gene Conversion #######
    ###############################   
    # Single Gene Conversion Events #
    ('blue', 'green', 'orange', 'orange') : 'N1G',
    ('blue', 'blue', 'orange', 'red') : 'N1L',  
    ('blue', 'orange', 'orange', 'purple') : 'N2G',    
    ('blue', 'blue', 'orange', 'yellow') : 'N2L',    
    ('blue', 'blue', 'orange', 'white') : 'N3G', 
    ('blue', 'empty', 'orange', 'orange') : 'N3L',
    
    # Double Gene Conversion Events #
    ('green', 'green', 'orange', 'orange') : 'O1G', 
    ('blue', 'blue', 'red', 'red') : 'O1L', 
        
    ('blue', 'orange', 'purple', 'white') : 'O2G', 
    ('blue', 'empty', 'orange', 'yellow') : 'O2L', 
        
    ('blue', 'green', 'orange', 'white') : 'O3G', 
    ('blue', 'empty', 'orange', 'red') : 'O3L',  
    
    ('blue', 'orange', 'orange', 'white') : 'O4G', 
        
    ('blue', 'blue', 'white', 'white') : 'O5G', 
    ('empty', 'empty', 'orange', 'orange') : 'O5L',  

    ('orange', 'orange', 'purple', 'purple') : 'O6G', 
    ('blue', 'blue', 'yellow', 'yellow') : 'O6L', 
    }


TETRAD_MAP_THREE = {
    ('blue','orange','orange'): 'A', 
    ('blue','blue','orange'): 'A', 
        
    ('blue','orange','purple'): 'B',
    ('orange','purple','yellow'): 'B', 
    ('blue','purple','yellow'): 'B',
    ('blue', 'orange' , 'yellow'): 'B', 
   
    ('blue','orange', 'white'): 'C',
    ('empty','orange', 'white'): 'C',  
    ('blue','empty', 'white'): 'C',
    ('blue','empty','orange'): 'C', 
        
    ('blue', 'orange', 'red'): 'D',
    ('blue', 'green', 'red'): 'D', 
    ('blue', 'green', 'orange'): 'D',  
    ('green','orange','red') : 'D',
        
    ('empty', 'orange', 'purple'): 'E', 
    ('empty','green', 'orange'): 'E', 
    ('empty','green', 'purple'): 'E', 
    ('green', 'orange', 'purple'): 'E',
        
    ('blue', 'white', 'yellow'): 'F', 
    ('blue', 'red', 'white'): 'F', 
    ('red', 'white', 'yellow'): 'F',  
    ('blue','red','yellow'): 'F', 
        
    ('purple', 'white', 'yellow'): 'G',
    ('empty', 'purple', 'white'): 'G',
    ('empty', 'purple', 'yellow'): 'G', 
    ('empty','white','yellow'): 'G', 
        
    ('purple', 'purple', 'yellow'): 'H', 
    ('purple', 'yellow', 'yellow'): 'H', 
        
    ('green', 'purple', 'red') : 'I',  
    ('green','purple','yellow') : 'I', 
    ('purple','red','yellow') : 'I', 
        
    ('empty','green','red') : 'J', 
    ('empty','green','white') : 'J', 
    ('green','red','white') : 'J',
    ('empty','red','white') : 'J', 
        
    ('green', 'red', 'red') : 'K', 
    ('green', 'green', 'red') : 'K',
        
    ('empty','empty', 'white'): 'L', 
    ('empty', 'white', 'white'): 'L',
        
    ('blue','blue', 'empty'): 'M',

    }


tetrad_counts = {}
def lookup_table(tetrad_type):
    tetrad_type = tuple(tetrad_type)
    
    # Update tetrad_counts
    if tetrad_type in tetrad_counts:
        tetrad_counts[tetrad_type] += 1
    else:
        tetrad_counts[tetrad_type] = 1

    # Check if tetrad is in the four-element mapping
    if tetrad_type in TETRAD_MAP_FOUR:
        return TETRAD_MAP_FOUR[tetrad_type], tetrad_counts 
    
    # Check if tetrad is in the three-element mapping
    three_element_tetrad = tetrad_type[:3]
    if three_element_tetrad in TETRAD_MAP_THREE:
        return TETRAD_MAP_THREE[three_element_tetrad], tetrad_counts 
    
    # Return default value if no match found
    return 'Z', tetrad_counts 


    

    
    
    
    
    
    
    
    
    
    
