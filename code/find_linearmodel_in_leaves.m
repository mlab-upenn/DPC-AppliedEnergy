function coeff = find_linearmodel_in_leaves(model, leafmodels, features_d)

NTrees = model.NTrees;
coeff = zeros(length(leafmodels{1}(1).coeff),NTrees);

for idt = 1:NTrees
    
    outputs = predict(model.Trees{idt}, features_d);

    % find the corresponding leaf index
    for idl = 1:length(leafmodels{idt})
        if abs(outputs - leafmodels{idt}(idl).mean)<1e-8
            break;
        end
    end

%     % litmus test
%     [testypred,testleaf] = predict(model.Trees{idt},features_d);
%     [~, testnode] = predict( model.Trees{idt}, model.X(~model.OOBIndices(:,idt),:) );
%     % if mean(model.Y(find(testnode==testleaf)))-testypred>1e-10
%     %     keyboard;
%     % end
%     if find(testnode==testleaf)~=leafmodels{idt}(idl).leaves{1,1}
%         keyboard;
%     end

    % extract the linear model coefficients from that leaf
    coeff(:,idt) = leafmodels{idt}(idl).coeff;
    
    if coeff(1,idt)<0 || coeff(1,idt)>50
        coeff(:,idt) = NaN;
    end
    
end
coeff(:,isnan(coeff(1,:))) = [];
coeff = mean(coeff,2);