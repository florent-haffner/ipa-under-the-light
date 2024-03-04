import torch
import torch.nn as nn


def create_dataloader(X, y, batch_size=16, shuffle=False):
    dataset = torch.utils.data.TensorDataset(X, y)
    return torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=shuffle)

def init_weights(m):
    if isinstance(m, nn.Conv1d):
        nn.init.kaiming_uniform_(m.weight, nonlinearity='leaky_relu')
    elif isinstance(m, nn.LazyLinear):
        nn.init.kaiming_uniform_(m.weight, nonlinearity='leaky_relu')
    elif isinstance(m, nn.MaxPool1d):
        nn.init.kaiming_uniform_(m.weight, nonlinearity='leaky_relu')
    return m

def kernel_regularization(model, loss, l2_lambda):
    """
    Regularization is only applied to Convolutions
    """
    l2_reg = 0
    # Simulate Keras's Kernel regularization
    layers_to_ignore = ['regressor']    
    for W in model.named_parameters():
         if "weight" in W[0]:
             layer_name = W[0].replace(".weight", "")
             # print(layer_name)
             if layer_name not in layers_to_ignore:
                 l2_reg = l2_reg + W[1].norm(2)

    loss = loss + l2_reg * l2_lambda
    return loss

"""
Train/test function
"""
def train(model, device, dataloader, l2_lambda, criterion, optimizer, lr_scheduler=None):
    model.train()  # Set the model to training mode
    last_loss = None
    for batch, (X, y) in enumerate(dataloader):
        optimizer.zero_grad()
        X, y = X.to(device), y.to(device)
        pred = model(X)
        loss = criterion(pred, y)

        loss = kernel_regularization(model, loss, l2_lambda)
        
        # Backward pass and optimization
        loss.backward()
        optimizer.step()
        
        if batch % 100 == 0:
            loss, current = loss.item(), (batch + 1) * len(X)
        
        last_loss = loss

    # Update the step of the scheduler
    if lr_scheduler: lr_scheduler.step()
    return last_loss.detach().cpu()

def test(model, device, dataloader, criterion):
    num_batches = len(dataloader)
    model.eval()
    test_loss = 0
    with torch.no_grad():
        for X, y in dataloader:
            X, y = X.to(device), y.to(device)
            pred = model(X)
            test_loss += criterion(pred, y).item()

    test_loss /= num_batches
    # return test_loss.detach().cpu()
    return test_loss
