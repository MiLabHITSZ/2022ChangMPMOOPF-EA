function WV = WV_generator2(N,m)
    H = 0;
    while nchoosek(H+m-1,m-1) <= N
        H = H+1;
    end
    H = H-1;
    n = nchoosek(H+m-1,m-1);
    
%     WV = zeros(n,m);
%     temp_WV =  nchoosek(1:H+m-1,m-1);
%     for i = 1:n
%         WV(i,:) = ([temp_WV(i,:)-[0:m-2]-1,H]-[0,temp_WV(i,:)-[0:m-2]-1])/H;    
%     end

%     WV = nchoosek(1:H+m-1,m-1) - repmat(0:m-2,n,1) - 1;
%     WV = ([WV,zeros(n,1)+H]-[zeros(n,1),WV])/H;

    WV_idx = nchoosek(1:(H+m-1),m-1);
    WV = zeros(n,m);
    for i = 1:n
        for j = 2:m-1
            WV(i,j) = (WV_idx(i,j)-WV_idx(i,j-1)-1)/H;
        end
    end
    WV(:,1) = (WV_idx(:,1)-1)/H;
    WV(:,m) = (H+m-1-WV_idx(:,m-1))/H;
    
    WV = max(WV,1e-6);
end