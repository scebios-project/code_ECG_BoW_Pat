function [feature_vote]=Vote(feature,Method,codebook,KNN,beta,lambda,sigma,A)
switch Method
    case {'VQ','SA','KCB','LLC','SPC'}
        feature_vote=zeros(size(feature,1),size(codebook,1));
        switch Method
            case {'VQ','SA','KCB'}
                disMatrix=dense_distance(codebook,feature,0,1,[],[],0,2);
                switch Method
                    case 'VQ'
                        [D,I]=min(disMatrix);
                        feature_vote(sub2ind(size(feature_vote),1:size(feature,1),I))=1;
                    case 'SA'
                        disMatrix=exp(-beta*disMatrix);
                        [D,I]=sort(disMatrix,'descend');
                        temp=bsxfun(@rdivide,D(1:KNN,:),sum(D(1:KNN,:),1));
                        for i=1:KNN
                        feature_vote(sub2ind(size(feature_vote),1:size(feature,1),I(i,:)))=temp(i,:);
                        end
                    case 'KCB'
                        disMatrix=(1/(sigma*(2*pi)^0.5))*exp(-0.5*(disMatrix.^2/(sigma^2)));
                        [D,I]=sort(disMatrix,'descend');
                        for i=1:KNN
                        feature_vote(sub2ind(size(feature_vote),1:size(feature,1),I(i,:)))=D(i,:);
                        end
                end
            case 'LLC'
                feature_vote=LLC_coding_appr(codebook,feature,KNN);
            case 'SPC'
%                 if isempty(A)
%                 beta = 1e-4;
%                 A = codebook*codebook' + 2*beta*eye(size(codebook,1));
%                 end
%                 Q = -codebook* feature';
%                 for i = 1:size(feature,1)
%                     feature_vote(i,:) = L1QP_FeatureSign_yang(lambda, A, Q(:,i));
%                 end;
                
                param.mode = 2;
                %param.mode = 1;
                param.K=size(codebook,1);
                param.L=size(feature,2);
                param.lambda=0.5;
                param.iter=100;
                %param.lambda=0.90;
                %param.iter=1000;
                %param.pos = true;
                
                feature_vote = mexLasso(feature',codebook',param);
                feature_vote = abs(feature_vote');
        end
    case 'FK'
        P=posterior(codebook,feature);
        feature_vote=zeros(size(feature,1),2*codebook.NDimensions*codebook.NComponents);
        for gmm_i=1:codebook.NComponents
            feature_vote(:,2*(gmm_i-1)*codebook.NDimensions+1:(2*gmm_i)*codebook.NDimensions)=...
                [(1/(size(feature,1)*codebook.PComponents(gmm_i)^0.5))*...
                bsxfun(@times,P(:,gmm_i),bsxfun(@rdivide,bsxfun(@plus,feature,codebook.mu(gmm_i,:)),codebook.Sigma(:,:,gmm_i))),...
                (1/(size(feature,1)*(2*codebook.PComponents(gmm_i))^0.5))*...
                bsxfun(@times,P(:,gmm_i),bsxfun(@rdivide,bsxfun(@plus,feature,codebook.mu(gmm_i,:)).^2,(codebook.Sigma(:,:,gmm_i).^2)-1))];
        end
    case 'SV'
end
end
