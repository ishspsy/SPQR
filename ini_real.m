function ini_b = ini_real(J,p,D,YY,tau)


ini_b2 = [];     %cell(1,J);  %zeros(2*J*p,1); %  cell(1,J);  % zeros(2*J*p,1);
    parfor kk = 1:J
    inds = union((2*kk-1):(2*J):(2*J*p),(2*kk):(2*J):(2*J*p));
    Xkk = D(kk:J:end,inds);
    Ykk = YY(kk:J:end);
    bkk =rq(Xkk, Ykk, tau);
    ini_b2 = [ini_b2; bkk];
    %ini_b2{kk} = bkk;
    %ini_b2(union((2*kk-1):(2*J):(2*J*p),(2*kk):(2*J):(2*J*p))) = bkk;
    end

 
 ini_b=[];
 parfor kk=1:J
 ini_b = [ini_b;reshape(ini_b2((2*(kk-1)*p+1):(2*kk*p)),2,p)];
 end
 