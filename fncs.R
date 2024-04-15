
expmod_win_prob <- function(k, l = c(2, 1, 3), a = c(.5, 1, 2)){
    
    # put constants in ascending order and rearrange other variables and recalculate index k
    
    m <- length(l)
    si <- sort(a, index.return=TRUE)
    
    k <- which(si$ix==k)
    a <- si$x
    d <- a-a[k]
    d <- c(d, Inf)
    l <- l[si$ix]
    
    I0 <- l[k] * exp(sum(l[1:(k - 1)]*d[1:(k - 1)]))/sum(l[1:k])*(1-exp(-sum(l[1:k])*d[k + 1]))
    
    I <- c()
    
    if(k == m) I <- 0
    
    else{
        for(j in 1:(m - k)){
            kk <- k + j
            ll <- l[1:kk][-k]
            dd<- d[1:kk][-k]
            I[j] <- l[k] * exp(sum(ll * dd)) /sum(l[1:kk])* (exp(-sum(l[1:kk])*d[kk]) - exp(-sum(l[1:kk])*d[kk + 1]))
            if(is.nan(I[j])) I[j] <- 0
        }
    }
    I0 + sum(I)
}

expmod_win_prob.vec <- Vectorize(expmod_win_prob, vectorize.args="k")



expmod_race_sim <- function(n=10000,l=c(1,2,3),a=c(.5,1,2)){
    zz <- matrix(0,nr=n,nc=length(l))
    for(i in 1:n){
        x <- rexp(length(l),l)
        y <- x+a
        zz[i,] <- sort(y, ind=T)$ix
    }
    zz
}




lratio <- function(l, i1 = NULL, i2) 
    ifelse(is.null(i1), 1/sum(l[i2]), sum(l[i1])/sum(l[i2]))


lexp <- function(l, a, i1 , i2) exp(- sum(l[i1])*(a[i2[1]] - a[i2[2]]))


r2 <- function(l, a){
    ifelse(a[1] > a[2], lratio(l, 1, 1:2) * lexp(l, a, 2, 1:2), 1 - lratio(l, 2, 1:2) * lexp(l, a, 1, 2:1))
}

r3 <- function(l, a){
    if(a[3] <= a[2]) {
        s1 <- r2(c(l[1], l[2]+l[3]), a[1:2]) * lratio(l, 2, 2:3) * lexp(l, a, 3,  2:3)
        s2 <- 0
    }
    
    else{
        s1 <- r2(c(l[1], l[2]+l[3]), a[c(1,3)]) * lratio(l, 2, 2:3) * lexp(l, a, 2,  3:2)
        
        if(a[3] < a[1]) s2 <- 0
        
        if((a[1] > a[2]) & (a[3] > a[1])) s2 <- lexp(l, a, 2, 1:2) * (lratio(l, 1, 1:2) * (1 - lexp(l, a, 1:2, c(3,1))) -
                                                                          (1 - lexp(l, a, 1, c(3, 1))) * lexp(l, a, 2, c(3,1)))
        
        
        if((a[2] > a[1]) & (a[3] > a[1])) s2 <- lexp(l, a, 1, 2:1) * (lratio(l, 1, 1:2) * (1 - lexp(l, a, 1:2, c(3,2))) -
                                                                          (1 - lexp(l, a, 1, c(3, 2))) * lexp(l, a, 2, c(3,2))) +
                (1 - lexp(l, a, 1, 2:1)) * (1 - lexp(l, a, 2, 3:2))
        
        
    }
    
    return(c(s1 = s1, s2 = s2, res = s1 + s2))
    
    
}


full_rank_probs_rec <- function(l, a){
    
    k <- length(l)
    
    if(k == 3) return(r3(l, a))
    
    if(a[k] <= a[k-1]) {
        new_l <- l[-k]
        new_l[k-1] <- l[k-1] + l[k]
        new_a <- a[-k]
        s1 <- full_rank_probs_rec(new_l, new_a)[3] * lratio(l, k-1, (k-1):k) * lexp(l, a, k , (k-1):k)
        s2 <- 0
    }
    
    else{
        new_l <- l[-k]
        new_l[k-1] <- l[k-1] + l[k]
        new_a <- a[-k]
        new_a[k-1] <- a[k]
        
        s1 <- full_rank_probs_rec(new_l, new_a)[3] * lratio(l, k-1, (k-1):k) * lexp(l, a, k-1,  k:(k-1))
        
        if(a[k] < a[k-2]) s2 <- 0
        
        if((a[k-2] > a[k-1]) & (a[k] > a[k-2])){
            new_l1 <- l[-k]
            new_l1[k-2] <- l[k-2] + l[k-1]
            new_l2 <- l[-k]
            s2 <- lexp(l, a, k-1, (k-2):(k-1)) * 
                (lratio(l, k-2, (k-2):(k-1)) * full_rank_probs_rec(new_l1, new_a)[2] - lexp(l, a, k-1, c(k, k-2)) *  
                     full_rank_probs_rec(new_l2, new_a)[2])
        }
        
        
        if((a[k-1] > a[k-2]) & (a[k] > a[k-1])) {
            new_l1 <- l[-k]
            new_l1[k-2] <- l[k-2] + l[k-1]
            new_l2 <- l[-k]
            new_a1 <- a[-k]
            new_a1[k-1] <- a[k]
            new_a1[k-2] <- a[k-1]
            new_a2 <- a[-k]
            s2 <- lexp(l, a, k-2, c(k-1, k-2)) *  
                (lratio(l, k-2, c(k-2, k-1)) * full_rank_probs_rec(new_l1, new_a1)[2] - lexp(l, a, k-1, c(k, k-1)) * 
                     full_rank_probs_rec(new_l2, new_a1)[2]) +
                full_rank_probs_rec(new_l2, new_a2)[2] * (1 - lexp(l, a, k-1, c(k, k-1)))
            
        }
        
    }
    return(c(s1 = s1, s2 = s2, res = s1 + s2))
}


full_rank_probs_rec_call <- function(l, a){
    full_rank_probs_rec(l, a)[3] %>% unname
}







rank_probs_rec_orig <- function(k, m=k, nracers = 5, rem = 1:nracers, l = runif(nracers), a = c(.2, .1, runif(3))){
    
    
    if(m == 2){
        
        res3 <- expmod_win2_prob(1:2, nracers, 1:nracers, l, a)
        
        return(res3)
        
    }
    
    
    else{
        
        if(a[1] > a[2]) {    
            
            l_new <- l[-1]
            a_new <- a[-1]
            a_new[1] <- a[1]
            a_new[2] <- a[3]
            
            l_new2 <- l[-1]
            l_new2[1] <- l[1] + l[2]
            
            res1 <- lexp(l, a, 2, 1:2) * (
                rank_probs_rec_orig(k, m - 1, nracers - 1, 1:(nracers - 1), l_new, a_new) - 
                    lratio(l, 2, 1:2) * rank_probs_rec_orig(k, m - 1, nracers - 1, 1:(nracers - 1), l_new2, a_new)
            )
            return(res1)
        }
        
        else{
            
            l_new <- l[-1]
            a_new <- a[-1]
            a_new[1] <- a[2]
            a_new[2] <- a[3]
            
            l_new2 <- l[-1]
            l_new2[1] <- l[1] + l[2]
            
            res2 <- rank_probs_rec_orig(k, m - 1, nracers - 1, 1:(nracers - 1), l_new, a_new) - 
                lexp(l, a, 1, 2:1) * lratio(l, 2, 1:2) * rank_probs_rec_orig(k, m - 1, nracers - 1, 1:(nracers - 1), l_new2, a_new)
            return(res2)
        }
    } 
}



rank_probs_rec <- function(k, m=k, nracers = 5, rem = 1:nracers, l = runif(nracers), a = c(.2, .1, runif(3))){
    
    
    if(m == 1){
        
        res <- expmod_win_prob(1, l, a)
        
        return(res)
        
    }
    
    else{
        
        if(a[1] > a[2]) {    
            
            l_new <- l[-1]
            a_new <- a[-1]
            a_new[1] <- a[1]
            a_new[2] <- a[3]
            
            l_new2 <- l[-1]
            l_new2[1] <- l[1] + l[2]
            
            res <- lexp(l, a, 2, 1:2) * (
                rank_probs_rec(k, m - 1, nracers - 1, 1:(nracers - 1), l_new, a_new) - 
                    lratio(l, 2, 1:2) * rank_probs_rec(k, m - 1, nracers - 1, 1:(nracers - 1), l_new2, a_new)
            )
            return(res)
        }
        
        else{
            
            l_new <- l[-1]
            a_new <- a[-1]
            a_new[1] <- a[2]
            a_new[2] <- a[3]
            
            l_new2 <- l[-1]
            l_new2[1] <- l[1] + l[2]
            
            res <- rank_probs_rec(k, m - 1, nracers - 1, 1:(nracers - 1), l_new, a_new) - 
                lexp(l, a, 1, 2:1) * lratio(l, 2, 1:2) * rank_probs_rec(k, m - 1, nracers - 1, 1:(nracers - 1), l_new2, a_new)
            return(res)
        }
    } 
}



expmod_win2_prob<- function(ws, nracers = 5, rem = 1:nracers, l = runif(nracers), a = runif(nracers)){
    
    if(a[ws[1]] < a[ws[2]]) {
        
        l_new <- l[-ws[1]]
        a_new <- a[-ws[1]]
        rem_new <- setdiff(rem, ws[1])
        ind_new <- which(rem_new == ws[2])
        
        s2 <- (1 - lexp(l, a, ws[1], rev(ws))) * expmod_win_prob(ind_new, l_new, a_new)
        
        l_new2 <- l_new
        l_new2[ind_new] <- sum(l[ws])
        
        l_new3 <- l[-ws[2]]
        rem_new2 <- setdiff(rem, ws[2])
        ind_new2 <- which(rem_new2 == ws[1])
        
        s1 <- lexp(l, a, ws[1], rev(ws)) * (
            expmod_win_prob(ind_new, l_new, a_new) - lratio(l, ws[2], ws) * expmod_win_prob(ind_new, l_new2, a_new)
        )
        
        res <- s1 + s2
    }
    
    else{
        
        l_new <- l[-ws[1]]
        a_new <- a[-ws[2]]
        rem_new <- setdiff(rem, ws[1])
        ind_new <- which(rem_new == ws[2])
        
        l_new2 <- l_new
        l_new2[ind_new] <- sum(l[ws])
        
        
        res <- lexp(l, a, ws[2], ws) * (expmod_win_prob(ind_new, l_new, a_new) -lratio(l, ws[2], ws) * expmod_win_prob(ind_new, l_new2, a_new))
    }
    return(res)
}


# =================================================================================================================================



gum_mod_probs  <- function(b = c(1, 2 ,3), a = c(0, 0, 0), alpha = 1){
    
    exp(-(b + a) / alpha)/sum(exp(-(b + a) / alpha))
    
}

race_nllh <- function(beta, model = model1, race_df, first_manche=TRUE, second_manche=TRUE, allranks = TRUE, maxrank = 10, nracers = 30){
    
    races <- unique(race_df$race)
    n_races <- length(races)
    
    if(model$type == "gumbel") eval(parse(text = model$alpha_ind))
    
    l1 <- 0
    
    if(first_manche){
        #first manche...
        
        for(race_number in 1:n_races){
            
            df <- subset(race_df, race == races[race_number])
            df <- arrange(df, time1) #arrange by first manche times
            nr <- nrow(df)
            
            
            eval(parse(text = model$lambda1))
            
            
            if(model$type == "gumbel")
            {
                for(i in 1:ifelse(allranks, nr - 1, maxrank)){
                    
                    l1 <- l1 - log(gum_mod_probs(lambda[i:nr], 0, alpha=exp(beta[alpha_ind]))[1])
                    
                }
            }
            
            else{
                
                m <- ifelse(allranks, nr, maxrank)
                l1 <- l1 - log(rank_probs_rec(m, m, nr, 1:nr, lambda, rep(0, nr)))
                
            }
        }
    }
    
    
    # total race...
    
    l2 <- 0
    
    if(second_manche){
        
        for(race_number in 1:n_races){
            df <- subset(race_df, race == races[race_number])
            df <- arrange(df, position) #arrange by final position
            nr <- nrow(df)
            
            eval(parse(text = model$lambda2))
            
            if(model$type == "gumbel")
            {
                for(i in 1:ifelse(allranks, nr - 1, maxrank)){
                    
                    l2<- l2 - log(gum_mod_probs(lambda[i:nr], df$time1[i:nr], alpha = exp(beta[alpha_ind]))[1])
                    
                }
            }
            
            else{
                
                m <- ifelse(allranks, nr, maxrank)
                l1 <- l1 - log(rank_probs_rec(m, m, nr, 1:nr, lambda, df$time1))
                
            }
        }
    }
    
    l <- l1 + l2 
    cat(l, fill=T)
    l
}    

fit <- function(model = model1, race_df = alldata, first_manche = TRUE, second_manche = TRUE, allranks = TRUE, maxrank = 10, nracers = 30, if_hessian = FALSE){

    n_races <- race_df$race %>% unique %>% length

    opt <- optim(rep(0, eval(parse(text=model$npar))), race_nllh, race_df = race_df, first_manche = first_manche, second_manche = second_manche, maxrank = maxrank, allranks = allranks, model = model, control = list(maxit = 5000),method="BFGS", hessian = if_hessian)
    
    if(if_hessian)
        se <- opt$hessian %>% solve %>% diag %>% sqrt
    else
        se <- NULL
    
    list(model = model, opt = opt, par = opt$par, se =se)
}



pred <- function(race_df, race_number, beta, model, plot = TRUE){
    
    races <- unique(race_df$race)
    n_races <- length(races)
    race_df <- arrange(race_df, bib2)
    
    df <- race_df[race_df$race == races[race_number],]
    nr <- nrow(df)
    
    eval(parse(text = model$lambda2))
    if(model$type == "gumbel") eval(parse(text = model$alpha_ind))
    
    if(model$type == "gumbel") 
        
        pr <-  gum_mod_probs(lambda, df$time1, alpha = exp(beta[alpha_ind]))
    
    else 
        
        pr <- expmod_win_prob.vec(1:nr, lambda, df$time1)
    
    p_df <- data.frame(surname=df$surname, bib=df$bib2, manche1_pos=df$rank1, points=df$points, prob=pr)
    
    if(plot == TRUE){
        p_plot <- ggplot(p_df, aes(bib,prob)) +geom_point(col='steelblue')+ xlab('Second Manche Start Position') + ylab('Win Probability') %>% print
        p_plot %>% print
    }
    
    p_df
}

#============================================================================================================================

# plotting functions

singleracerankplot <- function(race){
    data=read.csv(files[race])
    ggplot(data=data,aes(rank1,rank2)) + geom_point() + ggtitle('Position') + xlab("First Manche") + ylab("Second Manche") + 
        scale_x_continuous(breaks=c(1,5,10,15,20,25,30)) + scale_y_continuous(breaks=c(1,5,10,15,20,25,30)) + 
        stat_smooth(method='lm') +geom_abline(intercept=0, slope=1, colour="magenta")
}

labeller <- function(variable,value){
    return(races[value])
}



myColors <- brewer.pal(10, "Paired")
names(myColors) <- levels(alldata$race)
custom_colors <- scale_colour_manual(name = "Race", values = myColors)


rankplot <- function(data){
    ggplot(data=data,aes(rank1,rank2)) + geom_point(aes(colour=race)) + ggtitle('Position') + xlab("First Manche") + ylab("Second Manche") + 
        scale_x_continuous(breaks=c(1,5,10,15,20,25,30)) + scale_y_continuous(breaks=c(1,5,10,15,20,25,30)) + 
        stat_smooth(method='lm',colour="steelblue")+ labs(colour = "Event") + custom_colors
}

racerankplot <- function(data){
    ggplot(data=data,aes(rank1,rank2)) + geom_point(aes(colour=race)) + ggtitle('Position') + xlab("First Manche") + ylab("Second Manche") + 
        scale_x_continuous(breaks=c(1,5,10,15,20,25,30)) + scale_y_continuous(breaks=c(1,5,10,15,20,25,30)) + 
        stat_smooth(method='lm',colour="steelblue") + facet_wrap(race~., nrow=2) + guides(color = FALSE, size = FALSE)+ custom_colors
}



singlerankposplot <- function(data){
    ggplot(data=data,aes(rank1,position)) + geom_point(aes(colour=race)) + ggtitle('Position') + xlab("First Manche") + ylab("Overall") + 
        scale_x_continuous(breaks=c(1,5,10,15,20,25,30)) + scale_y_continuous(breaks=c(1,5,10,15,20,25,30)) + 
        stat_smooth(method='lm',colour="steelblue")+ labs(colour = "Event")+geom_abline(intercept=0, slope=1, colour="magenta")+ custom_colors
}

rankposplot <- function(data){
    ggplot(data=data,aes(rank1,position)) + geom_point(aes(colour=race)) + ggtitle('Position') + xlab("First Manche") + ylab("Overall") + 
        scale_x_continuous(breaks=c(1,5,10,15,20,25,30)) + scale_y_continuous(breaks=c(1,5,10,15,20,25,30)) + 
        stat_smooth(method='lm')+geom_abline(intercept=0, slope=1, colour="magenta")+ facet_wrap(race~., nrow=2) + guides(color = FALSE, size = FALSE)+ custom_colors
}

racetimeplot <- function(race){
    data=read.csv(files[race])
    ggplot(data=data,aes(time1,time2)) + geom_point() + ggtitle('Time (secs)') + xlab("First Manche") + ylab("Second Manche") +
        stat_smooth(method='lm')+ custom_colors
}

timeplot <- function(data){
    ggplot(data=data,aes(time1,time2)) + geom_point(aes(colour=race)) + ggtitle('Time (secs)') + xlab("First Manche") + ylab("Second Manche") +
        stat_smooth(method='lm',colour="steelblue") + labs(colour = "Event")  +geom_abline(intercept=0, slope=1, colour="magenta")+ custom_colors
}

septimeplot <- function(data){
    ggplot(data=data,aes(time1,time2)) + geom_point(aes(colour=race)) + ggtitle('Time (secs)') + xlab("First Manche") + ylab("Second Manche") +
        stat_smooth(method='lm',colour="steelblue")+ facet_wrap(race~., nrow=2, scales="free") + guides(color = FALSE, size = FALSE) + custom_colors
}

racetotaltimeplot <- function(race){
    data=read.csv(files[race])
    ggplot(data=data,aes(time1,time1+time2)) + geom_point() + ggtitle('Total time (secs)') + xlab("First Manche") + ylab("Total") +
        stat_smooth(method='lm',colour="steelblue")+ custom_colors
}

totaltimeplot <- function(data){
    ggplot(data=data,aes(time1,time1+time2)) + geom_point(aes(colour=race)) + ggtitle('Total time (secs)') + xlab("First Manche") + ylab("Total") +
        stat_smooth(method='lm',colour="steelblue")+ labs(colour = "Event")+ custom_colors
}

septotaltimeplot <- function(data){
    ggplot(data=data,aes(time1,time1+time2)) + geom_point(aes(colour=race)) + ggtitle('Total time (secs)') + xlab("First Manche") + ylab("Total") +
        stat_smooth(method='lm',colour="steelblue")+ facet_wrap(race~., nrow=2, scales="free")+ guides(color = FALSE, size = FALSE)+ custom_colors
}

