function v_ol = overlap(s1,s2,sigma,trSpan)
    v_ol = exp(-min([abs(s1-s2),trSpan-abs(s1-s2)])^2/(4*sigma^2));
end