function f = conformalmap(m,domain)

b = boundary(domain,'list');
bo = feval(m,b{1});
bi = feval(m,b{2});
f = conformalmap(@(z) feval(m,z),domain,region(bo,bi));
