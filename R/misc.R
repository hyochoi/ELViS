#' @import knitr
#' @noRd





candicol1 = c(brewer.pal(9,"Pastel1")[6], # candidate colors for exonic regions
              brewer.pal(8,"Pastel2")[6],
              brewer.pal(9,"YlOrBr")[1],
              brewer.pal(9,"YlOrRd")[1],
              brewer.pal(9,"YlOrRd")[2],
              brewer.pal(9,"Reds")[1],
              brewer.pal(9,"RdPu")[1],
              brewer.pal(9,"OrRd")[1],
              brewer.pal(9,"OrRd")[2],
              brewer.pal(9,"Oranges")[1],
              brewer.pal(9,"Oranges")[2],
              "aliceblue");
candicol2 = brewer.pal(12,"Set3") # candidiate colors for regions with shape changes
candicol3 = brewer.pal(8,"Pastel2") # candidiate colors for regions with shape changes
candicol = c(candicol2,candicol3);
exon.col = candicol1[9]



#' Get performance from the contingency table
#'
#' @param tbl 2 x 2 table with columns indicating predicted classes and rows indicating actual classes and the values indicating the number of cases. 1st row and column should indicate the target true event.
#'
#' @return performance metrics
#' @noRd
get_perf = function(tbl){
  rs = tbl |> apply(1,sum)
  cs = tbl |> apply(2,sum)
  Recall = Sensitivity = tbl[1,1]/rs[1]
  Specificity = tbl[2,2]/rs[2]
  Precision = PPV = tbl[1,1]/cs[1] # positive predictive value or precision
  NPV = tbl[2,2]/cs[2]
  F1 = (2*Precision*Recall)/(Precision+Recall)
  Accuracy = sum(diag(tbl))/sum(tbl)

  out = data.frame(
    Sensitivity = Sensitivity,
    Specificity = Specificity,
    PPV = PPV, # positive predictive value or precision
    NPV = NPV,
    Accuracy = Accuracy,
    F1 = F1
  )

  return(out)

}



#' Check if decimal places are trivially small, i.e. if it is integer
#'
#' @param x a numeric to test
#' @param th If `x - as.integer(x)` is less than th, this function returns TRUE
#'
#' @return logical indicating if input number is trivially small
#' @noRd
#' @examples
#'
#' is_zero_dec(1.3) # FALSE
#' is_zero_dec(1.01) # FALSE
#' is_zero_dec(1.3-1) # TRUE
#'
is_zero_dec = function(x,th=10^(-10)){
  (x - as.integer(x)) < th
}




#' Internal function for get_A_B_Prop_Subcount
#'
#' @param A A count
#' @param A_str A name
#' @param B B count
#' @param B_str B name
#' @param target_prop desired proportions of A (= A/(A + B))
#' @param n_cycle_th maximum number of revision cycles
#'
#' @return list containint updated A and B that satisfies target proportion
#' @noRd
get_A_B_Prop_Subcount__internal =
  function(A,A_str,B,B_str,target_prop,n_cycle_th = 10){
    n_cycle = 1
    while( n_cycle <= n_cycle_th ){

      print(paste0("n_cycle : ",n_cycle))

      A = B*target_prop/(1-target_prop)
      if(is_zero_dec(A)){break}
      print(A_str)
      A = as.integer(A)

      B = A/target_prop - A
      if(is_zero_dec(B)){break}
      print(B_str)
      B = as.integer(B)

      n_cycle = n_cycle + 1
    }
    return(
      structure(list(A,B),
                names = c(A_str,B_str)
                )

      )
  }


#' Get an exact integer number pair with starting numbers A and B that satisfies target proportion
#'
#' @param A target number A
#' @param B target number B
#' @param target_prop target proportion the user want to achieve in the form of A/(A+B).
#' @param n_cycle_th number of adjustment cycle to make both A and B integer-valued
#'
#' @return list containint updated A and B that satisfies target proportion
#' @noRd
#'
#' @examples
#' N_alt_ori = 50L
#' N_ctrl_ori = 600L
#' get_A_B_Prop_Subcount(A=N_alt_ori,B=N_ctrl_ori,target_prop=0.01)
#' get_A_B_Prop_Subcount(A=N_alt_ori,B=N_ctrl_ori,target_prop=0.05)
#' get_A_B_Prop_Subcount(A=N_alt_ori,B=N_ctrl_ori,target_prop=0.1)
#' get_A_B_Prop_Subcount(A=N_alt_ori,B=N_ctrl_ori,target_prop=0.2)
get_A_B_Prop_Subcount =
  function(A,B,target_prop,n_cycle_th = 10){
    A_str = deparse(substitute(A))
    B_str = deparse(substitute(B))

    ori_prop = A/(A + B)
    if(ori_prop > target_prop){
      output =
        get_A_B_Prop_Subcount__internal(A=A,A_str = A_str,B=B,B_str = B_str,target_prop = target_prop,n_cycle_th=n_cycle_th)
    }else if(ori_prop < target_prop){
      output =
        get_A_B_Prop_Subcount__internal(A=B,A_str = B_str,B=A,B_str = A_str,target_prop = 1-target_prop,n_cycle_th=n_cycle_th) %>%
        rev
    }else{
      print("No change")
      output = structure(list(A,B),
                         names = c(A_str,B_str))
    }
    return(output)
  }
