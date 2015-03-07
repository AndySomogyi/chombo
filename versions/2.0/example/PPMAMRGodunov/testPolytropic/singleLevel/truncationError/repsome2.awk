#
# Get the "nvar"th component name
#
# The input is assumed to contain lines of the form:
#
#   component-name: L1-norm, L2-norm, Linf-norm
#
# for each component.  It can contain other lines not of this form.
#

# Start counting component names
BEGIN {
  cvar = 0;
}

# Look for lines with the above format - there will be one for each component
NF == 4 && $1 ~ ".*:" && $2 ~ ".*," && $3 ~ ".*," {
  # The current component number
  cvar++;

  # If this equals the component number asked for then output its name
  if (cvar == nvar) {
    name = substr($1,1,index($1,":")-1);
    printf("%s",name);
    exit;
  }
}
