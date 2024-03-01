TETRAD_MAP_FOUR = {
    ## Non- Crossover ##
    ('blue','blue','orange','orange'): '1',
    ##########################
    ### Single Crossovers ####
    ##########################
    ('blue','orange','purple','yellow'): '2',  
    ('blue','empty','orange', 'white'): '3',
    ###########################
    #### Double Crossovers ####
    ###########################
    ('blue','green', 'orange', 'red'): '4', 
    ('empty','green', 'orange', 'purple'): '5',
    ('blue','red', 'white', 'yellow'): '6',
    ('empty','purple', 'white', 'yellow'): '7',
    ('purple','purple', 'yellow', 'yellow'): '8',    
    ###########################
    #### Triple Crossovers ####
    ###########################
    ('green','purple','red','yellow') : '9',
    ('empty','green','red','white') : '10',
    ###########################
    #### Quadra Crossovers ####
    ###########################
    ('green','green', 'red', 'red') : '11',
    ###########################
    ######### MI and MII NDJ ##########
    ###########################
    ('empty','empty', 'white', 'white'): '12',
    ('blue','blue', 'empty', 'orange') : '13', 
    ('blue','empty','orange','orange') : '14',
    ###############################
    ######  Gene Conversion #######
    ###############################   
    # Single Gene Conversion Events #
    ('blue', 'green', 'orange', 'orange') : 'N1G',
    ('blue', 'blue', 'orange', 'red') : 'N1L',  
    ('blue', 'orange', 'orange', 'purple') : 'N2G',    
    ('blue', 'blue', 'orange', 'yellow') : 'N2L',    
    ('blue', 'blue', 'orange', 'white') : 'N3G', 
    
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
    }


TETRAD_MAP_THREE = {
    ('blue','orange','orange'): '1', 
    ('blue','blue','orange'): '1', 
        
    ('blue','orange','purple'): '2',
    ('orange','purple','yellow'): '2', 
    ('blue','purple','yellow'): '2',
    ('blue', 'orange' , 'yellow'): '2', 
   
    ('blue','orange', 'white'): '3',
    ('empty','orange', 'white'): '3',  
    ('blue','empty', 'white'): '3',
    ('blue','empty','orange'): '3', 
        
    ('blue', 'orange', 'red'): '4',
    ('blue', 'green', 'red'): '4', 
    ('blue', 'green', 'orange'): '4',  
    ('green','orange','red') : '4',
        
    ('empty', 'orange', 'purple'): '5', 
    ('empty','green', 'orange'): '5', 
    ('empty','green', 'purple'): '5', 
    ('green', 'orange', 'purple'): '5',
        
    ('blue', 'white', 'yellow'): '6', 
    ('blue', 'red', 'white'): '6', 
    ('red', 'white', 'yellow'): '6',  
    ('blue','red','yellow'): '6', 
        
    ('purple', 'white', 'yellow'): '7',
    ('empty', 'purple', 'white'): '7',
    ('empty', 'purple', 'yellow'): '7', 
    ('empty','white','yellow'): '7', 
        
    ('purple', 'purple', 'yellow'): '8', 
    ('purple', 'yellow', 'yellow'): '8', 
        
    ('green', 'purple', 'red') : '9',  
    ('green','purple','yellow') : '9', 
    ('purple','red','yellow') : '9', 
        
    ('empty','green','red') : '10', 
    ('empty','green','white') : '10', 
    ('green','red','white') : '10',
    ('empty','red','white') : '10', 
        
    ('green', 'red', 'red') : '11', 
    ('green', 'green', 'red') : '11',
        
    ('empty','empty', 'white'): '12', 
    ('empty', 'white', 'white'): '12',
        
    ('blue','blue', 'empty'): '13',

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


    

    
    
    
    
    
    
    
    
    
    
