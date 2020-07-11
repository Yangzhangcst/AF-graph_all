function [Z] = unifiedcluster(K,alpha,beta)
K = K*K';
maximum = max(K(:));
K = K ./maximum;

addpath('KSC/qpc')
[~,n]=size(K);
Z=eye(n);
for i=1:200
    Zold=Z;
    Z= (Z+Z')/2;
    D = diag(sum(Z));
    L = D-Z; 
    F=eig1(L, n, 0);

    parfor ij=1:n
        all = veccomp2(ij,n,F);
        all=all.*all;
        H=2*alpha*eye(n)+2*K;
        H=(H+H')/2;
        ff=beta/2*all'-2*K(:,ij);
        % we use the free package to solve quadratic equation: http://sigpromu.org/quadprog/index.html
        Z(:,ij) = qpas(H,ff,[],[],ones(1,n),1,zeros(n,1),ones(n,1));
    end
    err = norm(Z-Zold)/norm(Zold);
    if err<1e-5,  break;  end

end
end

function [all]=veccomp2(ij,n,F)
    for ji=1:n
        all(ji)=norm(F(ij,:)-F(ji,:));
    end
end
