if(!require(pomp)){install.packages('pomp');library(pomp)}

write.sweep = function(
  reps=10,
  epsilon_neg = c(0,1),
  epsilon_pos = c(0,1),
  lambda = c(.1,.9),
  mu = c(-5,-2),
  founders = c(1,25),
  max_founder_date = c(7,365*3),
  R0.local = c(1,1),
  SUP_rate = c(-2,0),
  SUP_max = c(2,2),
  numLoci = c(10,25),
  MOI.source = c(1,4),
  geno.error = c(0,.1)){
  scalars = sobol(vars = list(
    epsilon_neg = epsilon_neg,
    epsilon_pos = epsilon_pos,
    lambda = lambda,
    mu = mu,
    founders = founders,
    max_founder_date = max_founder_date,
    R0.local = R0.local,
    # SUP_rate = SUP_rate,
    # SUP_max = SUP_max,
    numLoci = numLoci,
    MOI.source = MOI.source,
    geno.error = geno.error),
    reps)
  for(ii in 1:reps){
    write.sim.data(
      file = ii,
      epsilon_neg = scalars$epsilon_neg[ii],
      epsilon_pos = scalars$epsilon_pos[ii],
      lambda = scalars$lambda[ii],
      mu = 10 ^ scalars$mu[ii],
      founders = round(scalars$founders[ii]),
      max_founder_date = round(scalars$max_founder_date[ii]),
      R0.local = scalars$R0.local[ii],
      SUP_toggle = FALSE,
      # SUP_rate = 10 ^ scalars$SUP_rate[ii],
      # SUP_max = round(scalars$SUP_max[ii]),
      max_cases = 50,
      numLoci = round(scalars$numLoci[ii]),
      MOI.source = scalars$MOI.source[ii],
      geno.error = scalars$geno.error[ii])
  }
}


write.sim.data = function(
  file = 1,
  epsilon_neg = 1e-6,
  epsilon_pos = 1e-6,
  lambda = .5,
  mu = 1e-6,
  founders = 35,
  max_founder_date = 500,
  R0.local = 1,
  SUP_toggle = FALSE,
  SUP_rate = .5,
  SUP_max = 2,
  max_cases = 75,
  max_tries = 20,
  shape.serial = 22.09,
  scale.serial = 2.19,
  numLoci = 10,
  MOI.source = 3,
  geno.error = .1,
  # For negative binomial offspring distribution
  k = 2)
{
  # define function to check for cycles
  if(!require(igraph)){install.packages('igraph'); library(igraph)}
  
  repeat{
    # simulate network
    edges = matrix(1,founders,2)
    edges[,2] = (1:founders)+1
    dates.infection = runif(founders) * max_founder_date
    
    tovisit = (1:founders)+1
    
    while(length(tovisit) > 0){
      parent = tovisit[1]
      #number.offspring = rpois(1,R0.local)
      number.offspring = rnbinom(n=1,size=k,prob=k/(R0.local+k))
      
      if(number.offspring > 0){
        tovisit = c(tovisit,(max(edges)+1):(max(edges)+number.offspring))
        edges = rbind(
          edges,
          cbind(
            rep(parent,number.offspring),
            (max(edges)+1):(max(edges)+number.offspring)))
        dates.infection = c(dates.infection, dates.infection[parent-1] + pmax(rgamma(number.offspring,shape=shape.serial,scale=scale.serial),rep(0,number.offspring)))
      }
      
      tovisit = tovisit[-1]
      
      if(max(edges) >= max_cases + 1)
        break
      
      print("in tovisit")
    }
    print("out tovisit")
    
    network = NULL
    network = matrix(0,nrow=max(edges),ncol=max(edges))
    network[edges] = 1
    
    if(nrow(network) == max_cases + 1){break}
  }
  
  
  # write simulated network to file
  string.to.write = ''
  network.ind = which(network > 0, arr.ind = T)
  for(ii in 1:length(which(network > 0))){
    if(network.ind[ii,1] == 1){
      string.to.write = paste(
        string.to.write,
        's-',
        as.character(network[network.ind[ii,1],network.ind[ii,2]]), '-',
        as.character(network.ind[ii,2]-1), ';', sep = '')
    }
    if(network.ind[ii,1] > 1){
      string.to.write = paste(
        string.to.write,
        as.character(network.ind[ii,1]-1), '-',
        as.character(network[network.ind[ii,1],network.ind[ii,2]]), '-',
        as.character(network.ind[ii,2]-1), ';', sep = '')
    }
  }
  write.table(
    string.to.write,
    file=paste('../data/tmp/network_true_',file,'.csv',sep=''),
    quote=FALSE,sep='',row.names=FALSE,col.names=FALSE)
  
  
  # fit parameters to simulate allele frequency distributions
  f = list(
    Ara2 = as.numeric(read.table('../data/prevalences_empirical/frequency_Ara2.txt',header=FALSE,sep='\t',quote='')[-1]),
    PFG377 = as.numeric(read.table('../data/prevalences_empirical/frequency_PFG377.txt',header=FALSE,sep='\t',quote='')[-1]),
    PfPK2 = as.numeric(read.table('../data/prevalences_empirical/frequency_PfPK2.txt',header=FALSE,sep='\t',quote='')[-1]),
    PolyAlpha = as.numeric(read.table('../data/prevalences_empirical/frequency_PolyAlpha.txt',header=FALSE,sep='\t',quote='')[-1]),
    TA1 = as.numeric(read.table('../data/prevalences_empirical/frequency_TA1.txt',header=FALSE,sep='\t',quote='')[-1]),
    TA40 = as.numeric(read.table('../data/prevalences_empirical/frequency_TA40.txt',header=FALSE,sep='\t',quote='')[-1]),
    TA60 = as.numeric(read.table('../data/prevalences_empirical/frequency_TA60.txt',header=FALSE,sep='\t',quote='')[-1]),
    TA81 = as.numeric(read.table('../data/prevalences_empirical/frequency_TA81.txt',header=FALSE,sep='\t',quote='')[-1]),
    TA87 = as.numeric(read.table('../data/prevalences_empirical/frequency_TA87.txt',header=FALSE,sep='\t',quote='')[-1]))
  max_freqs = sapply(f,function(ff)max(ff))
  min_freqs = sapply(f,function(ff)min(ff))
  library(MASS)
  min_freq_allowable = 0.05
  f.fit = fitdistr(
    (unlist(f)[unlist(f)>=min_freq_allowable]-min_freq_allowable)/(1-min_freq_allowable),
    'beta', start=list(shape1=0.5, shape2=0.5))$estimate
  
  # simulated genetics are defined by four features
  numAlleles = rpois(numLoci, mean(c(11,8,17,21,26,21,12,13,14)))
  freqs.source = list()
  max_sim_freqs = numeric()
  min_sim_freqs = numeric()
  for(i in 1:numLoci){
    freqs.source[[i]] = rbeta(numAlleles[i],f.fit[1],f.fit[2])*(1-min_freq_allowable)+min_freq_allowable
    max_sim_freqs[i] = max(freqs.source[[i]])
    min_sim_freqs[i] = min(freqs.source[[i]])
  }
  
  # declare object to store simulated alleles for each individual
  alleles = list()
  
  # assign each founder an MOI
  moi = rpois(founders, MOI.source)
  while(sum(moi == 0) | sum(moi > 6)){moi = rpois(founders,MOI.source)}
  
  # loop over each locus ll
  for(ll in 1 : length(freqs.source)){
    
    # assign a temporary matrix with a row for each sample and column for each allele
    alleles.temp = matrix(0, nrow = nrow(network) - 1, ncol = length(freqs.source[[ll]]))
    
    # loop over each founder
    for(ff in which(network[1,] > 0)){
      # draw alleles for ff according to their frequencies in the source population
      # alleles.temp[ff-1,sample(ncol(alleles.temp),moi[ff-1],replace=T,prob=freqs.source[[ll]])] = 1
      for(aa in 1:ncol(alleles.temp)){
        alleles.temp[ff-1,aa] = rbinom(1,1,freqs.source[[ll]][aa])
      }
    }
    
    # loop over each allele aa at the locus
    # draw alleles along each edge of the network
    for(ee in (founders+1):nrow(edges)){
      locus.save = alleles.temp[edges[ee,2]-1,]
      locus.temp = rep(0,length(freqs.source[[ll]]))
      while(sum(locus.temp) == 0){
        locus.temp = rep(0,length(freqs.source[[ll]]))
        for(aa in 1 : length(freqs.source[[ll]])){
          if(alleles.temp[edges[ee,1]-1, aa] > 0){
            locus.temp[aa] = rbinom(1, 1, lambda)
          }else if(alleles.temp[edges[ee,1]-1, aa] == 0){
            locus.temp[aa] = rbinom(1, 1, mu)
          }  
        }
        print("in")
        print(alleles.temp[edges[ee,1]-1,])
        print(lambda)
      }
      print("out")
  
      alleles.temp[edges[ee,2]-1,] = locus.save + locus.temp
      rm(locus.temp)
    }
    
    alleles.temp[alleles.temp > 1] = 1
    
    alleles[[ll]] = alleles.temp  
  }

  # simulate false positives and false negatives  
  prop_eps_pos = runif(length(alleles), min=0, max=1)
  prop_eps_neg = runif(length(alleles), min=0, max=1)

  for(i in 1:length(alleles)){
    for(nn in 1:nrow(alleles[[i]])){
      alleles.save = alleles[[i]][nn,]
      repeat{
        alleles[[i]][nn,] = alleles.save
        for(aa in 1:length(alleles[[i]][nn,])){
          if(alleles[[i]][nn,aa] > 0){
            allele_freq = freqs.source[[i]][aa]
            eps_neg_allele = prop_eps_neg * (1-allele_freq)
            if(rbinom(1,1,eps_neg_allele)){
              alleles[[i]][nn,aa] = 0
            }              
          }else if(alleles[[i]][nn,aa] == 0){
            allele_freq = freqs.source[[i]][aa]
            eps_pos_allele = prop_eps_pos * allele_freq
            if(rbinom(1,1,eps_pos_allele)){
              alleles[[i]][nn,aa] = 1
            }
          }
        }   
        if(sum(alleles[[i]][nn,]) > 0){
          break
        }
      }
    }  
  }
  
  # simulate locus-wide genotyping failure
  for(ll in 1:length(alleles)){
    for(nn in 1:nrow(alleles[[ll]])){
      if(runif(1) < geno.error){
        alleles[[ll]][nn,] = rep(0,ncol(alleles[[ll]]))
      }
    }
  }
  
  
  # write true parameter values to file
  write.table(
    paste(
      paste('epsilon_neg_L', 1:numLoci, ',', prop_eps_neg, sep='', collapse='\n'),
      '\n',
      paste('epsilon_pos_L', 1:numLoci, ',', prop_eps_pos, sep='', collapse='\n'),
      '\n',
      'lambda,', lambda, '\n',
      'mu,', mu, '\n',
      'founders,', founders, '\n',
      'max_founder_date,', max_founder_date, '\n',
      'R0.local,', R0.local, '\n',
      'numLoci,', numLoci, '\n',
      'MOI.source,', MOI.source, '\n',
      'geno.error,', geno.error, sep=''),
    file=paste('../data/tmp/scalars_true_',file,'.csv',sep=''),
    quote=FALSE,sep='',row.names=FALSE,col.names=FALSE)
  
  # write data about nodes to file
  write(
    paste(
      'Time',
      paste('L',1:length(alleles),sep='',collapse=','),sep=','),
    paste('../data/tmp/nodes_',file,'.csv',sep=''))
  for(nn in 1:(nrow(network)-1)){
    write(
      paste(
        floor(dates.infection[nn]),
        paste(sapply(alleles,function(aa)
          paste(aa[nn,],collapse='')),collapse=','),sep=','),
      paste('../data/tmp/nodes_',file,'.csv',sep=''),
      append=TRUE)
  }
  
  
  # write frequencies to file
  write(
    paste(
      's',
      paste('L',1,sep=''),
      paste(freqs.source[[1]],collapse=','),sep=','),
    paste('../data/tmp/sources_',file,'.csv',sep=''))
  for(ll in 2:length(freqs.source)){
    write(
      paste(
        's',
        paste('L',ll,sep=''),
        paste(freqs.source[[ll]],collapse=','),sep=','),
      paste('../data/tmp/sources_',file,'.csv',sep=''),
      append=TRUE)
  }
  for(ll in 1:length(freqs.source)){
    write(
      paste(
        'local',
        paste('L',ll,sep=''),
        paste(pmax(colMeans(alleles[[ll]]),rep(.00001,ncol(alleles[[ll]]))),collapse=','),sep=','),
      paste('../data/tmp/sources_',file,'.csv',sep=''),
      append=TRUE)
  }
}


isCyclic = function(edges)
{
  network = matrix(0,nrow=max(edges),ncol=max(edges))
  network[edges] = 1
  g = graph.adjacency(
    network,
    mode = 'directed')
  
  return(!is.dag(g))
}
