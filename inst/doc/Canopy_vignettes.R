### R code from vignette source 'Canopy_vignettes.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Installation1 (eval = FALSE)
###################################################
## install.packages('Canopy')


###################################################
### code chunk number 2: Installation2 (eval = FALSE)
###################################################
## install.packages("devtools")
## library(devtools)
## install_github("yuchaojiang/Canopy/package")


###################################################
### code chunk number 3: Input
###################################################
library(Canopy)
data("MDA231")

projectname = MDA231$projectname ## name of project
R = MDA231$R; R ## mutant allele read depth (for SNAs)
X = MDA231$X; X ## total depth (for SNAs)
WM = MDA231$WM; WM ## observed major copy number (for CNA regions)
Wm = MDA231$Wm; Wm ## observed minor copy number (for CNA regions)
epsilonM = MDA231$epsilonM ## standard deviation of WM, pre-fixed here
epsilonm = MDA231$epsilonm ## standard deviation of Wm, pre-fixed here
## Matrix C specifices whether CNA regions harbor specific CNAs 
## only needed if overlapping CNAs are observed, specifying which CNAs overlap
C = MDA231$C; C
Y = MDA231$Y; Y ## whether SNAs are affected by CNAs


###################################################
### code chunk number 4: Tree_elements1
###################################################
data('MDA231_tree')
MDA231_tree$Z # Z matrix specifies the position of the SNAs along the tree branch
MDA231_tree$cna.copy # major and minor copy number (interger values) for each CNA
MDA231_tree$CM # Major copy per clone for each CNA
MDA231_tree$Cm # Minor copy per clone for each CNA
MDA231_tree$Q # whether an SNA precedes a CNA


###################################################
### code chunk number 5: Tree_elements2
###################################################
MDA231_tree$H # if an SNA precedes a CNA, whether it resides in the major copy
MDA231_tree$P # clonal compostion for each sample
MDA231_tree$VAF # VAF based on current tree structure


###################################################
### code chunk number 6: Sampling1 (eval = FALSE)
###################################################
## K = 3:6 # number of subclones
## numchain = 20 # number of chains with random initiations
## sampchain = canopy.sample(R = R, X = X, WM = WM, Wm = Wm, epsilonM = epsilonM, 
##             epsilonm = epsilonm, C = C, Y = Y, K = K, numchain = numchain, 
##             simrun = 50000, writeskip = 200, projectname = projectname,
##             cell.line = TRUE, plot.likelihood = TRUE)
## save.image(file = paste(projectname, '_postmcmc_image.rda',sep=''),
##            compress = 'xz')


###################################################
### code chunk number 7: Sampling2
###################################################
data("MDA231_sampchain")
sampchain = MDA231_sampchain
k = 3
K = 3:6
sampchaink = MDA231_sampchain[[which(K==k)]]


###################################################
### code chunk number 8: Sampling3
###################################################
length(sampchain) ## number of subtree spaces (K=3:6)
length(sampchain[[which(K==4)]]) ## number of chains for subtree space with 4 subclones
length(sampchain[[which(K==4)]][[1]]) ## number of posterior trees in each chain


###################################################
### code chunk number 9: BIC
###################################################
burnin = 100
thin = 10
bic = canopy.BIC(sampchain = sampchain, projectname = projectname, K = K,
               numchain = numchain, burnin = burnin, thin = thin, pdf = FALSE)
optK = K[which.max(bic)]


###################################################
### code chunk number 10: fig1
###################################################
# Note: this segment is soley for generating BIC plot in the vignettes.
# Use Canopy.BIC() with pdf = TRUE to generate this plot directly.
par(mfrow=c(1,2))
projectname = 'MDA231'
numchain = 20
clikelihood = matrix(nrow = numchain, ncol = length(sampchaink[[1]]), data = NA)
for(numi in 1:numchain){
  for(i in 1:ncol(clikelihood)){
    clikelihood[numi,i] = sampchaink[[numi]][[i]]$likelihood
  }
}
plot(1:ncol(clikelihood), clikelihood[1,], type='l', xlab = 'Iteration',
     ylab = 'Log-likelihood', col = 1, ylim = c(min(clikelihood), 
                                                max(clikelihood)))
for(numi in 2:numchain){
  points(1:ncol(clikelihood), clikelihood[numi,], type = 'l', col = numi)
}
title(main=paste('Posterior likelihood', k, 'clones', numchain,
            'chains'),cex=0.6)
plot(K, bic, xlab = 'Number of subclones', ylab = 'BIC', type = 'b', xaxt = "n")
axis(1, at = K)
abline(v = (3:6)[which.max(bic)], lty = 2)
title('BIC for model selection')


###################################################
### code chunk number 11: Post
###################################################
post = canopy.post(sampchain = sampchain, projectname = projectname, K = K,
                 numchain = numchain, burnin = burnin, thin = thin, optK = optK,
                 C = C, post.config.cutoff = 0.05)
samptreethin = post[[1]]   # list of all post-burnin and thinning trees
samptreethin.lik = post[[2]]   # likelihoods of trees in samptree
config = post[[3]] # configuration for each posterior tree
config.summary = post[[4]] # configuration summary
print(config.summary)
# first column: tree configuration
# second column: posterior configuration probability in the entire tree space
# third column: posterior configuration likelihood in the subtree space


###################################################
### code chunk number 12: Plot
###################################################
config.i = config.summary[which.max(config.summary[,3]),1]
cat('Configuration', config.i, 'has the highest posterior likelihood!\n')
# plot the most likely tree in the posterior tree space
output.tree = canopy.output(post, config.i, C)
canopy.plottree(output.tree)

# plot the tree with configuration 1 in the posterior tree space
output.tree = canopy.output(post, 1, C)
canopy.plottree(output.tree,pdf=TRUE,pdf.name = 
                    paste(projectname,'_first_config.pdf',sep=''))


###################################################
### code chunk number 13: fig2
###################################################
canopy.plottree(output.tree)


###################################################
### code chunk number 14: Try it your self (eval = FALSE)
###################################################
## library(Canopy)
## data(toy)
## projectname = 'toy'
## R = toy$R; X = toy$X; WM = toy$WM; Wm = toy$Wm
## epsilonM = toy$epsilonM; epsilonm = toy$epsilonm; Y = toy$Y
## 
## K = 3:6; numchain = 10
## sampchain = canopy.sample(R = R, X = X, WM = WM, Wm = Wm, epsilonM = epsilonM, 
##                           epsilonm = epsilonm, C = NULL, Y = Y, K = K, 
##                           numchain = numchain, simrun = 50000, writeskip = 200,
##                           projectname = projectname, cell.line = FALSE,
##                           plot.likelihood = TRUE)


###################################################
### code chunk number 15: fig3
###################################################
data(toy)
canopy.plottree(toy$besttree, txt = FALSE, pdf = FALSE)


###################################################
### code chunk number 16: Try it your self2 (eval = FALSE)
###################################################
## library(Canopy)
## data(toy2)
## projectname = 'toy2'
## R = toy2$R; X = toy2$X; WM = toy2$WM; Wm = toy2$Wm
## epsilonM = toy2$epsilonM; epsilonm = toy2$epsilonm; Y = toy2$Y
## true.tree = toy2$true.tree
## 
## K = 3:6; numchain = 10
## sampchain = canopy.sample(R = R, X = X, WM = WM, Wm = Wm, epsilonM = epsilonM, 
##                           epsilonm = epsilonm, C = NULL, Y = Y, K = K, 
##                           numchain = numchain, simrun = 50000, writeskip = 200,
##                           projectname = projectname, cell.line = FALSE,
##                           plot.likelihood = TRUE)


###################################################
### code chunk number 17: fig4
###################################################
data(toy2)
canopy.plottree(toy2$true.tree, txt = FALSE, pdf = FALSE)


###################################################
### code chunk number 18: sessionInfo
###################################################
toLatex(sessionInfo())


