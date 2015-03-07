#
# Get the "norm" error norm for component "compname"
#
# The input is assumed to contain lines of the form:
#
#   component-name: L1-norm, L2-norm, Linf-norm
#
# for each component.  It can contain other lines not of this form.
#

# Look for lines with the above format - there will be one for each component
NF == 4 && $1 ~ ".*:" && $2 ~ ".*," && $3 ~ ".*," {
  curname = substr($1,1,index($1,":")-1);

  if (curname == compname) {
    if (norm == "L1") {
      error = substr($2,1,index($2,",")-1);
      printf("%17.10e",error);
    } else
    if (norm == "L2") {
      error = substr($3,1,index($3,",")-1);
      printf("%17.10e",error);
    } else
    if (norm == "Linf") {
      error = $4;
      printf("%17.10e",error);
    } else {
      print "Huh?"
    }
  }
}
