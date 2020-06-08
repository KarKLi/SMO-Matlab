syms t;
alpha=[1+0*t,0*t,0*t,0*t,1+0*t];
y=[1,1,1,-1,-1];
x=[1,2;2,3;3,3;2,1;3,2];
if sum(alpha.*y)~=0
    fprintf("总和为%.3f\n",sum(alpha.*y));
    fprintf("总和不为0，数据错误\n");
else
    % SMO 算法
    % SMO algorithm
    % 首先，我们尝试性地给数据+t
    % Firstly, we add variable t to alpha.
    W_value = W(alpha,y,x); 
    fprintf("W的最大值为：%.4f\n",W_value);
    fprintf("\n");
    for i=1:size(y,2)
        for j=1:size(y,2)
            if (j == i)
                continue
            end
            alpha_new = alpha;
            if y(i)*y(j) > 0 % 代表i,j同号，一个+t一个-t
                % It means i and j have the same positive/negative symbol
                % We assume that +t to alpha(i) and -t to alpha(j)
                alpha_new(i) = alpha_new(i)+t;
                alpha_new(j) = alpha_new(j)-t;
                if sum(alpha.*y)~=0
                    fprintf("总和为%.3f\n",sum(alpha.*y));
                    fprintf("总和不为0，数据错误\n");
                    quit(1);
                end
                % 找出最好的t
                % Solving the best t
                t_best = solve(diff(W(alpha_new,y,x)==0),t);
                % 此时要检验t_best会不会令alpha有小于0的分量，若有，缩放该t
                % Before adding the t, we need to check if there is an alpha which
                % will become zero after change.
                if alpha(i) + t_best < 0 || alpha(j)-t_best < 0
                    % 这个时候要看t_best的正负，若t_best为正，那么就要缩放成alpha(j)
                    % Checking the pos/nega of t_best
                    % if t_best is positive, t_best = alpha(j)
                    if t_best > 0
                        t_best = alpha(j);
                    else
                        % 如果t_best为负，就要缩放成alpha(i)
                        % if t_best is negative, t_best = alpha(i)
                        t_best = alpha(i);
                    end
                end
                alpha(i) = alpha(i)+t_best;
                alpha(j) = alpha(j)-t_best;
            else % 代表i,j异号，两个都+t
                % It means i and j have the same positive/negative symbol
                % +t to both alpha(i) and alpha(j)
                alpha_new(i) = alpha_new(i)+t;
                alpha_new(j) = alpha_new(j)+t;
                if sum(alpha.*y)~=0
                    fprintf("总和为%.3f\n",sum(alpha.*y));
                    fprintf("总和不为0，数据错误\n");
                    quit(1);
                end
                % 找出最好的t
                % Solving the best t
                t_best = solve(diff(W(alpha_new,y,x)==0),t);
                % 此时要检验t_best会不会令alpha有小于0的分量，若有，缩放该t
                % Before adding the t, we need to check if there is an alpha which
                % will become zero after change.
                if alpha(i) + t_best < 0 || alpha(j)+t_best < 0
                    % 这个时候t_best要被缩放到令一个量刚好为0
                    % 假如t_best = -4,alpha(i) = 3,alpha(j) = 2
                    % 那t_best应该要被缩放成-2
                    % t_best needs to be changed to some alpha equals to zero.
                    t_best = max(-alpha(i),-alpha(j));
                end
                alpha(i) = alpha(i)+t_best;
                alpha(j) = alpha(j)+t_best;
            end
            % alpha此时已被更新，计算W
            % alpha has be updated, calculating the new W(alpha)
            alpha
            W_value = W(alpha,y,x);
            fprintf("W的最大值为：%.4f\n",W_value);
            fprintf("\n");
        end
    end
    fprintf("W的最大值为：%.4f\n",W_value);
    % 回求最优w*
    % w是一个长度为size(x,2)的行向量
    % Calculating the w* with formula
    % \sum{alpha_i*y_i*x_i}
    w = zeros(1,size(x,2));
    for i=1:size(y,2)
        w = w + alpha(i)*y(i)*x(i,:);
    end
    % 求b*
    % 找到负例中w*^Txi最大的值
    % 在这里为了看着方便，w是行向量，而x也是行向量，所以与原公式有点出入，但结果没错
    % Calculating the b*
    % b* satisfies this formula
    % -(max_{i:y^{(i)}=-1}w^*^Tx^{(i)}+min_{i:y^{(i)}=1}w^*^Tx^{(i)})/2
    % Be cautious: For personal habitation, w is row vector and so does x
    % So the real calculation will be a litter different from original formula.
    max_value = 0;
    min_value = 1e8;
    for i=1:size(y,2)
        if y(i) == -1
            max_value = max(max_value,w*x(i,:)');
        else
            min_value = min(min_value,w*x(i,:)');
        end
    end
    b = -(max_value+min_value)/2;
    fprintf("最优的w为：");
    w
    fprintf("最优的b为：");
    b
end

function [res] = W(alpha,y,x)
    res = sum(alpha);
    for i=1:5
        for j = 1:5
            res = res - 0.5*(alpha(i)*alpha(j)*y(i)*y(j)*x(i,:)*x(j,:)');
        end
    end
end
