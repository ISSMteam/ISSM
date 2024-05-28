% function to replace Nan by 0 in J, then rescale with the mass matrix
%
function J = rescalegradientNan(md, unscaledJ)

nanflag = isnan(unscaledJ);
unscaledJ(nanflag) = 0;
J = rescalegradient(md, unscaledJ);
