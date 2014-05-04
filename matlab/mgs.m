% append column onto QR decomp with modified Gram-Schmidt
function [A,Q,R,zt,rhs] = mgs(A,Q,R,t,zt,rhs)

   if t == 1
      Q(:,t) = A(:,t);
      R(t,t) = norm(Q(:,t),2);

      if abs(R(t,t)) < 10^(-14)
         error('dun goofed')
      end
      Q(:,t) = Q(:,t) / R(t,t);

      zt = dot(Q(:,t),rhs);
      rhs = rhs - zt*Q(:,t);
   
   else
      
      for i=1:t-1
         R(i,t) = dot(Q(:,i),A(:,t));
         A(:,t) = A(:,t) - R(i,t)*Q(:,i);
      end

      Q(:,t) = A(:,t);
      R(t,t) = norm(Q(:,t),2);

      if abs(R(t,t)) < 10^(-14)
         error('dun goofed')
      end
      Q(:,t) = Q(:,t) / R(t,t);

      zt = dot(Q(:,t),rhs);
      rhs = rhs - zt*Q(:,t);
 
   end

end
