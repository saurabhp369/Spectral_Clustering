data = load('data1.mat');
X = data.X;
M = 8;
symbols = ['x','o','+', '_', '*', '^', 'd', 's'];
clusters = {X};
while(length(clusters) < M)
    new_cluster = {};
    for c = 1:length(clusters)
        N = length(clusters{:,c});
        W = zeros(N);
        for i = 1:N
            for j = 1:N
                if i == j
                    continue
                end
                W(i,j) = exp(-norm(clusters{:,c}(:,i) - clusters{:,c}(:,j))^2/10);
            end
        end
        
        D = zeros(N);
        for i = 1:N
            D(i,i) = sum(W(:,i));
        end
        
        [v, w] = eigs(D^-0.5 *(D-W)*D^-0.5,2,'smallestabs');
        z2 = v(:,2);
        
        labels = sign(D^-0.5 * z2);
        idx1 = find(labels == -1);
        idx2 = find(labels == 1);
        new_cluster(end+1) = {clusters{:,c}(:,idx1)};
        new_cluster(end+1) = {clusters{:,c}(:,idx2)}; 
    end
    clusters = new_cluster;
    l = 0;
    for b = 1:length(clusters)
        loss = 0;
        points = clusters{:,b};
        x = points(1,:);
        y = points(2,:);
        for u = 1:length(points)
            for v = 1:length(points)
                loss = loss + norm(points(:,u)-points(:,v))^2;
            end
        end
        l = l + loss/length(points);
        scatter(x,y,symbols(b))
        axis off
        hold on
    end
    J = 0.5*l
    chr = int2str(length(clusters));
    filename = append('spectral data1 ', chr);
    t = append('classes-',chr )
    title(t)
    saveas(gcf,filename,'png')
    clf  
end