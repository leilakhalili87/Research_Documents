function [mo,count] = equi_metric(n)
tol1 = 1e-12;
count = 0;
for a=-n:n
    for b=-n:n
        for c=-n:n
            for d=-n:n
                for e=-n:n
                    for f=-n:n
                        for g=-n:n
                            for h=-n:n
                                if (a*e ~= b*d)
                                    ip = ( 1-(b*f*g + c*d*h - c*e*g -a*f*h))/(a*e - b*d);
                                    if (abs(round(ip) - ip) < tol1)
                                        if (abs(ip) <= n+tol1)
                                            count=count+1;
                                            w=[[a,b,c];[d,e,f];[g,h,ip]];
                                            mo{count}=w;
                                        end
                                    end
                                    ip = (-1-(b*f*g + c*d*h - c*e*g -a*f*h))/(a*e - b*d);
                                    if (abs(round(ip) - ip) < tol1)
                                        if (abs(ip) <= n+tol1)
                                            count=count+1;
                                            w=[[a,b,c];[d,e,f];[g,h,ip]];
                                            mo{count}=w;
                                        end
                                    end
                                else
                                    if (abs(b*f*g + c*d*h - c*e*g -a*f*h)-1 == 0)
                                        for ip=-n:n
                                            count=count+1;
                                            w=[[a,b,c];[d,e,f];[g,h,ip]];
                                            mo{count}=w;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
end