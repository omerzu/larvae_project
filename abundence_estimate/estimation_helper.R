###################
#Optimization######
###################

###################
#Helper Functions#
###################

#@in X,Y N*O matrices each row represent one sample, 
#X = the reads, Y - areas
#@out vector w (w*O) of the estimated wjs
remove.degeneracy = function(values) {
    return(values[which(values != 0 & values !=1)]);
}


learnW = function (X,Y) {
    Xs = colSums(X);
    Ys = colSums(Y);
    W = Ys/Xs;
    W[which(is.nan(W))] = 1;
    return(W) 
}
learnW = function (X,Y,lambda=0) {
    W = rep(0,dim(X)[2]);
    N = dim(X)[1];
    for (i in 1:dim (X)[2]) {
        xi = X[,i];
        yi = Y[,i];
        W[i] = (lambda + sum(yi*xi))/(lambda + sum(xi*xi));
    }
    W[which(is.nan(W))] = 1;
    return(W) 
}
predict = function (Xi,W){
    return(Xi*W / sum(Xi*W));
}
predictReads = function (testData, W, catg.factors){
    indices = sapply(testData$catg,function(x,catg){which(x == catg)},catg.factors);
    nreads = testData$Reads * W[indices] / sum(testData$Reads * W[indices]);
    #nreads = testData$Reads * W[indices];
    return(nreads);
}
predictReads.round = function (testData,W,catg.factors){
    indices = sapply(testData$catg,function(x,catg){which(x == catg)},catg.factors);
    nreads = round(testData$Reads * W[indices])
    return(nreads);
}
predictReads.order = function (testData,W,catg.factors){
    indices = sapply(testData$Order,function(x,catg){which(x == catg)},catg.factors);
    nreads = testData$Reads * W[indices] / sum(testData$Reads * W[indices]);
    return(nreads);
}
predictReads.order.round = function (testData,W,catg.factors){
    indices = sapply(testData$Order,function(x,catg){which(x == catg)},catg.factors);
    nreads = round(testData$Reads * W[indices]);
    return(nreads);
}

###List functions###
sort.assignment.list=function(kBestAssignments){
  list.len = length(kBestAssignments) / 2;
	for (i in 1:list.len) {
		for (j in (i+1):list.len) {
			if (j>list.len) {next;}
			elm_i_score= kBestAssignments[[2*i]];
			elm_j_score= kBestAssignments[[2*j]];
			if(elm_i_score < elm_j_score) {#Swap elements
				tmp_elm = kBestAssignments[[2*i-1]];
				kBestAssignments[[2*i-1]] = kBestAssignments[[2*j-1]];
				kBestAssignments[[2*i]] = elm_j_score;
				kBestAssignments[[2*j-1]] = tmp_elm;
				kBestAssignments[[2*j]] = elm_i_score
			}
		}
	}
  return(kBestAssignments);
}
updateKassingments = function (kBestAssignments,k,cur_assignemt,score) {
  for(i in 1:k){ #check if this assignment exists
    if(sum(cur_assignemt==kBestAssignments[[2*i-1]])==length(cur_assignemt)) {
      return(kBestAssignments);
    }
  }
  kBestAssignments[[2*k-1]] = cur_assignemt;
	kBestAssignments[[2*k]] = score;
	kBestAssignments = sort.assignment.list(kBestAssignments);
  return(kBestAssignments);
}
###get the population vector and returens s*2 array of species and it's comulative area##
calculateArea = function (pop,s) {
  #pop - a f*2 array, the first column is the fish sp. id and the second is it's area
  #s in the #species
  areas = cbind(1:s,rep(0,length(s)))
  for(i in 1:(dim(pop)[1])) {
    areas[pop[i,1],2] = areas[pop[i,1],2] + pop[i,2];
  }
  return(areas);
}
###Takes an assignment and returns its log-likelihood (without constants)##
calculateAssignmentLikelihood = function (assignment,reads,s) {
  #assignment - f*2 vector, each row contain fish id and it's area
  #reads s*1 vector, each row contains the reads fraction of sp.i
  # - s #species
  areas = calculateArea(assignment,s);
  areas[,2] = areas[,2] / sum(areas[,2]);
  llhood = 0;
  llhood = sum(reads*log(areas[,2]))
  return(llhood)
}
calculateAssignmentLikelihood.L2norm = function (assignment,reads,s) {
  #assignment - f*2 vector, each row contain fish id and it's area
  #reads s*1 vector, each row contains the reads fraction of sp.i
  # - s #species
  areas = calculateArea(assignment,s);
  areas[,2] = areas[,2] / sum(areas[,2]);
  readsFrac = reads/sum(reads)
  score = norm(readsFrac - areas[,2],type="2");
  return(-score)
}
###Given assignment sample a new assignment which is likely better###
neighbourAssignmnet = function(assignment,readsFraction,s) {
  f = dim(assignment)[1]
  cur_areas = calculateArea(pop=assignment,s=s);
  cur_areas = cur_areas[,2] / sum (cur_areas[,2]);
  expected_areas = readsFraction;
  change_probs  = rep (1,s);
  pick_probs  = rep (1,s);
  #calculte fish picking probabilits
  change_probs = change_probs + (cur_areas - expected_areas); #increse sampling prob for over estimeted populations
  pick_probs = pick_probs - (cur_areas - expected_areas);
  #notmalize the probs:
  change_probs = change_probs / sum(change_probs);
  pick_probs = pick_probs / sum(pick_probs);
  #for each fish: decide if to change assignment w.r.t probs 
  isChanged = F;
  while(!isChanged){
    for (i in 1:f) {
      if (runif(n=1) < change_probs[assignment[i,1]]) { #change assigment
        assignment[i,1] = sample(x=1:s,size=1,prob=pick_probs);
        isChanged = TRUE;
      }
   }
  }
  return(assignment)
}
randomAssignmnet = function(readsFraction,s,f,isUniform=FALSE) {
  assignment = array(rep(NA,2*f),dim=c(f,2));
  #pick one fish for each species
  assignment[sample(1:f,size=s,replace=FALSE),1] = 1:s
  indexToFill = which(is.na(assignment[,1]));
  assignment[indexToFill,1] = sample(1:s,size = length(indexToFill),prob=readsFraction,replace=TRUE)
  if(isUniform) {
    assignment[indexToFill,1] = sample(1:s,size = length(indexToFill),replace=TRUE);
  }
  return(assignment);
}
###Manage the optimizatoin process###
#Assignemnt that doesn't assign any fish to one of the species are invalid. In those cases we need to pick larvae to missing species
makeAssignmentValid = function (assignment,s) {
    if (length(unique(assignment[,1])) == s) {
      return(assignment); #assignemnt is OK
    }
    else { #make sure we picked all the species
      f = dim(assignment)[1];
      assignment[sample(1:f,size=s,replace=FALSE),1] = 1:s;
      return(assignment);
    }
}
###Pick new assignment based on the current one###
neighbourAssignmnet = function(assignment,reads,s) {
  f = dim(assignment)[1]
  readsFraction = reads / sum(reads);
  cur_areas = calculateArea(pop=assignment,s=s);
  cur_areas = cur_areas[,2] / sum (cur_areas[,2]);
  expected_areas = readsFraction;
  change_probs  = rep (1,s);
  pick_probs  = rep (1,s);
  #calculte fish picking probabilits
  change_probs = change_probs + (cur_areas - expected_areas); #increse sampling prob for over estimeted populations
  pick_probs = pick_probs - (cur_areas - expected_areas);
  #notmalize the probs:
  change_probs = change_probs / sum(change_probs);
  pick_probs = pick_probs / sum(pick_probs);
  #for each fish: decide if to change assignment w.r.t probs 
  isChanged = F;
  while(!isChanged){
    for (i in 1:f) {
      if (runif(n=1) < change_probs[assignment[i,1]]) { #change assigment
        assignment[i,1] = sample(x=1:s,size=1,prob=pick_probs) ;
        #change_probs[assignment[i,1]] = change_probs[assignment[i,1]] / 2; #reduce the chance to resample this type
        #change_probs = change_probs / sum(change_probs); #re-normalize;
        isChanged = TRUE;
      }
   }
  }
  return(makeAssignmentValid(assignment,s))
}

searchBestAssignemt = function(reads,fish_areas,max_trials = 1000000,seed =0) {
#observedFractionOfReads - s*1 vector, the ith element is the reads fraction of sp. i
#fish_areas - f*1 vector, the ith element contains the area of the ith fish
set.seed(seed);
observedFractionOfReads  =reads / sum(reads);
s = length(observedFractionOfReads);
f = length(fish_areas);
IDs = sample(x=1:s,size=f,replace=TRUE)
best_assignment = cbind(IDs,fish_areas);
best_assignment = makeAssignmentValid(best_assignment,s);
cur_assignemt = best_assignment;
best_llhood = calculateAssignmentLikelihood (best_assignment,reads,s);
llhood = best_llhood;
assign_cnt = 0;
steps_cnt = 0;
while (assign_cnt < max_trials) {
	cur_assignemt = neighbourAssignmnet(best_assignment,observedFractionOfReads,s);
	llhood = calculateAssignmentLikelihood (cur_assignemt,reads,s);
  if (llhood > best_llhood) {
      best_assignment = cur_assignemt;
      best_llhood = llhood;
	  assign_cnt = round(assign_cnt * 0.90); #get extra moves for each improvment
    }
  assign_cnt = assign_cnt + 1;
  steps_cnt = steps_cnt + 1;
}
  return(list(best_assignment,best_llhood,seed,steps_cnt));
}
par.searchBestAssignemt = function (reads,fish_areas,max_trials,initalSeed,seedsNumber,isL2=FALSE) {
  seeds = ((initalSeed):(seedsNumber+initalSeed));
  if(isL2) {
    results = lapply(X=seeds,FUN=searchBestAssignemt.L2norm,reads=reads,fish_areas=fish_areas,max_trials=max_trials);
  }
  else{
  results = lapply(X=seeds,FUN=searchBestAssignemt,reads=reads,fish_areas=fish_areas,max_trials=max_trials);
  }
  best_asn = results[[1]][[1]];
  best_lhood = results[[1]][[2]];
  best_seed = results[[1]][[3]]
  best_steps = results[[1]][[4]]
  for (i in 1:seedsNumber) {
    r = results[[i]];
    if(r[[2]] > best_lhood) {
      best_asn = r[[1]]; 
      best_lhood = r[[2]];
      best_seed = r[[3]];
	  best_steps = r[[4]];
    }
  }
  return(list(best_asn,best_lhood,best_seed,best_steps))
}
searchBestAssignemt.L2norm = function(reads,fish_areas,max_trials=1000000,seed=0 ) {
#observedFractionOfReads - s*1 vector, the ith element is the reads fraction of sp. i
#fish_areas - f*1 vector, the ith element contains the area of the ith fish
set.seed(seed);
observedFractionOfReads  =reads / sum(reads);
s = length(observedFractionOfReads);
f = length(fish_areas);
IDs = sample(x=1:s,size=f,replace=TRUE)
best_assignment = cbind(IDs,fish_areas);
best_assignment = makeAssignmentValid(best_assignment,s);
cur_assignemt = best_assignment;
best_llhood = calculateAssignmentLikelihood.L2norm (best_assignment,reads,s);
llhood = best_llhood;
sec_best_assigment = best_assignment;
sec_best_llhood = best_llhood;
assign_cnt = 0;
while (assign_cnt < max_trials) {
  llhood_exp =  exp(best_llhood - llhood);
  if(is.na(llhood_exp)){ llhood_exp = 0};
  cur_assignemt = neighbourAssignmnet(best_assignment,observedFractionOfReads,s);
  llhood = calculateAssignmentLikelihood.L2norm(cur_assignemt,reads,s);
  if (llhood > best_llhood) {
    if (!(sum(best_assignment[,1] == cur_assignemt[,1])==(dim(best_assignment)[2]))) { #make sure it's not the same assignment that was sampled
        sec_best_assigment = best_assignment;
      }
      sec_best_assigment = best_assignment
      best_assignment = cur_assignemt;
      sec_best_llhood = best_llhood;
      best_llhood = llhood;
    }
  assign_cnt = assign_cnt + 1;
}
  return(list(best_assignment,best_llhood,seed));
}
#Currently the scaling is already done in imageJ
searchKBestAssignemt.L2norm = function(reads,fish_areas,max_trials = 1000000,k=1) {
#observedFractionOfReads - s*1 vector, the ith element is the reads fraction of sp. i
#fish_areas - f*1 vector, the ith element contains the area of the ith fish
observedFractionOfReads  =reads / sum(reads);
s = length(observedFractionOfReads);
f = length(fish_areas);
IDs = sample(x=1:s,size=f,replace=TRUE)
assignment_score=NA;
best_assignment = cbind(IDs,fish_areas);
best_assignment = makeAssignmentValid(best_assignment,s);
cur_assignemt = best_assignment;
best_llhood = calculateAssignmentLikelihood.L2norm (best_assignment,reads,s);
llhood = best_llhood;
prev_assignment = cur_assignemt;
kBestAssignments = list();
for(i in 1:k) {
 kBestAssignments[[2*i-1]]=best_assignment;
 kBestAssignments[[2*i]]=best_llhood;
}
minKscore = best_llhood
assign_cnt = 0;
while (assign_cnt < max_trials) {
  llhood_exp =  exp(best_llhood - llhood);
  if(is.na(llhood_exp)){ llhood_exp = 0};
  cur_assignemt = neighbourAssignmnet(best_assignment,observedFractionOfReads,s);
  llhood = calculateAssignmentLikelihood.L2norm(cur_assignemt,reads,s);
  if (llhood > best_llhood) {
    if (!(sum(best_assignment[,1] == cur_assignemt[,1])==(dim(best_assignment)[2]))) { #make sure it's not the same assignment that was sampled
        sec_best_assigment = best_assignment;
      }
      sec_best_assigment = best_assignment
      best_assignment = cur_assignemt;
      sec_best_llhood = best_llhood;
      best_llhood = llhood;
    }
  if(sum(cur_assignemt != prev_assignment)){
    if(llhood  > minKscore) {
      kBestAssignments = updateKassingments(kBestAssignments,k,cur_assignemt,llhood);
      minKscore = kBestAssignments[[2*k]];
    }
  }
  assign_cnt = assign_cnt + 1;
  prev_assignment=cur_assignemt;
}
  return(list(best_assignment,best_llhood,kBestAssignments));
}

##simulation of different reads number
simulateMultinomialReads = function(readsFactor=1) {
  data  =readDataFile();
  data$catg=data$Family;#no normalization, work by family
  catg.factors = unique(x=data$catg);
  samples = unique(data$Sample);
  repeats = 50;
  res.colnames = c("sample","identifedFraction","numOfFamilies","numOfLarvae");
  results = matrix(0,nrow=repeats*length(samples),ncol=length(res.colnames));
  colnames(results) = res.colnames;
  for (i in 1:length(samples)) {
    sampleNum = samples[i];
    sampleData = data[which(data$Sample == sampleNum),];
    lData = convertRawDataToTrainingData(sampleData,catg.factors)
    families.areas = lData[[1]];
    sample.reads = lData[[2]];
    sample.areas = sampleData$Area;
    sample.areas = sample.areas / sum(sample.areas);
    totalReads=sum(sampleData$Reads)* readsFactor;
    numOfLarvae = length(sample.areas);
    s = length(remove.degeneracy(families.areas));
    test.famlies = unique(as.character(sampleData$Family));
    test.true.assignment = cbind(as.numeric(sapply(X=sampleData$Family,FUN=function(x,families){which(x==families)},test.famlies)),sample.areas);
    #print (paste("sample=",sampleNum,"families=",s,"numOfLarvae=",numOfLarvae));
    for (r in 1:repeats) {
      sampledReads = c();
      while(length (sampledReads) != s){ #make sure to sample all species
      sampledReads = sample(x=1:s,size=totalReads,replace=TRUE,prob=remove.degeneracy(families.areas));
      sampledReads = table(sampledReads);
      sampledReads = sampledReads / sum(sampledReads);
      }
      test.trails = 10000;
      lAssignmentResults = searchBestAssignemt.L2norm(sampledReads,sample.areas,test.trails);
      best_asn = lAssignmentResults[[1]];
      results[r+(i-1)*repeats,1] = sampleNum;
      results[r+(i-1)*repeats,2] = scoreAssignment.fishwise(best_asn,test.true.assignment) / numOfLarvae;
      results[r+(i-1)*repeats,3] = s;
      results[r+(i-1)*repeats,4] = numOfLarvae;
    }
  }
  return(as.data.frame(results));
}


