function y = trim_silence(x)
    env = abs(x); thresh = 0.03*max(env);
    idx = find(env>thresh);
    if isempty(idx), y=x(:); return; end
    margin = round(0.1*length(x));
    y = x(max(1,idx(1)-margin) : min(length(x),idx(end)+margin));
end
