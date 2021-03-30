# Function to compute the CNE estimate of the total number of species
# We assume that x is a species by samples matrix of the number of specimens
# of each species in each sample
cneFxn <- function( comm )
{
    # Safeguard if only one species or one sample 
    if ( dim(comm)[1] == 0 || dim(comm)[2] == 0 )
    {
        print( "comm is not a matrix" )
        print( comm )
        return( nrow( comm ) )
    }

    # Number of species in each sample
    numSpecies <- colSums( comm > 0 )

    numSingletons <- colSums( comm == 1 )
    numDoubletons <- colSums( comm == 2 )

    y  <- comm > 0
    y1 <- y[ rowSums(y) == 1, ]
    y2 <- y[ rowSums(y) == 2, ]

    # We need to make special provisions for y1 or y2 not being tables.
    # If y1 or y2 is one row, we simply multiply by 1 to convert logical to int, and we are done.
    # If y1 or y2 is empty, then we need a vector of same length as there are rows
    if ( ncol(y) == 1 )                     # Single site, all are singletons
        numSingletonSpecies <- nrow(y)
    else if ( length( y1 ) > ncol(y) )      # We have a table, typical case
        numSingletonSpecies <- colSums( y1 )
    else if ( length( y1 ) == ncol(y) )     # Single species being a singleton
        numSingletonSpecies <- y1 * 1
    else                                    # No species being singletons
        numSingletonSpecies <- rep(0, times=ncol(y))

    if ( ncol(y) == 1 )                     # Single site, no doubletons
        numDoubletonSpecies <- 0
    else if ( length( y2 ) > ncol(y) )      # We have a table, typical case
        numDoubletonSpecies <- colSums( y2 )
    else if ( length( y1 ) == ncol(y) )     # Single species being a doubleton
        numDoubletonSpecies <- y1 * 1
    else                                    # No species being doubletons
        numDoubletonSpecies <- rep(0, times=ncol(y))
 
    estimatedSingletonSpecies <- numSingletonSpecies
    estimatedDoubletonSpecies <- numDoubletonSpecies

    for ( i in 1:length(estimatedSingletonSpecies) )
    {
        estimatedSingletonSpecies[i] <- estimatedSingletonSpecies[i] + 
            ((numSingletons[i]*(numSingletons[i]-1))/(2.0*(numDoubletons[i]+1)))*(numSingletonSpecies[i]/numSpecies[i])
    } 
 
    for ( i in 1:length(estimatedDoubletonSpecies) )
    {
        estimatedDoubletonSpecies[i] <- estimatedDoubletonSpecies[i] + 
            ((numSingletons[i]*(numSingletons[i]-1))/(2.0*(numDoubletons[i]+1)))*(numDoubletonSpecies[i]/numSpecies[i])
    } 
    
    Q1 <- sum( estimatedSingletonSpecies )
    Q2 <- 0.5 * sum( estimatedDoubletonSpecies )
    
    nrow(comm) + Q1*(Q1-1) / (2.0 * (Q2 + 1))
}

