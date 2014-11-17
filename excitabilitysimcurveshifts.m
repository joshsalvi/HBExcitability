%load('');


% Shift the curves in order to calculate the difference between them at a
% later time.
tr2=find(pkspikeratedet(1,:)==0);
tr2=max(tr2);
for j = 1:maxiter
    pkspikshift2(j) = (transitionend(j)) - (tr2);
    pkspikeratestoAVG_shift(j,:) = pkspikeratestoAVG(j,:);
    pkspikeratestoAVG_shift(j,1:length(pkspikeratestoAVG(j,abs(pkspikshift2(j)):end))) = pkspikeratestoAVG(j,abs(pkspikshift2(j)):end);
    pkspikeratestoAVG_shift(j,length(pkspikeratestoAVG(j,abs(pkspikshift2(j)):end)):end) = 0;
end

% Find the difference between curves
for j = 1:maxiter
    for k = 1:maxiter
        % Find the difference between shifted curves (curves)
        pkspikeratestoAVG_diff{j,k} = pkspikeratestoAVG(j,:) - pkspikeratestoAVG(k,:);
        pkspikeratestoAVG_diffshift{j,k} = pkspikeratestoAVG_shift(j,:) - pkspikeratestoAVG_shift(k,:);
        % Normalized difference between shifted curves (one value)
        % [diff/sum]
        pkspikeratestonormsum(j,k) = (sum(pkspikeratestoAVG(j,1:tr2+pkspikshift2(j)))-sum(pkspikeratestoAVG(k,1:tr2+pkspikshift2(k))))/(sum(pkspikeratestoAVG(j,1:tr2+pkspikshift2(j)))+sum(pkspikeratestoAVG(k,1:tr2+pkspikshift2(k))));
    end
end

% Find widths of transition zones at 1-99% of maximum
for j = 1:maxiter
for k = 1:99
    maxj = pkspikeratestoAVG(j,(pkspikeratestoAVG(j,1:transitionend(j)) == max(pkspikeratestoAVG(j,1:transitionend(j)))));
    tr3 = findnearest(pkspikeratestoAVG(j,:),max1*0.01*k);tr3=tr3(1);
    transitionwidth(j,k) = I(transitionend(j))-I(tr3);
    if k == 50
        transitionmid(j) = I(tr3);  % Control parameter at which spike rate is 50% of maximum
    end
end
end
    
    
