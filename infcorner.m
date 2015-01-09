function zeta = infcorner(angle1,angle2)

zeta = [ homog(exp(1i*angle1),0), homog(exp(1i*angle2),0) ];

end