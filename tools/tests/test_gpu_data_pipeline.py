import pytest
import torch
from copy import deepcopy
from torch.utils.data import DataLoader, TensorDataset

from scripts.runtime import (
    assert_nested_equal,
    move_to_device,
    optimized_loader_options,
    scalar_history_to_list,
    seeded_generator,
)


def fixed_dataset():
    features = torch.arange(64 * 8, dtype=torch.float32).reshape(64, 8) / 100
    labels = torch.arange(64, dtype=torch.long) % 2
    return TensorDataset(features, labels)


def test_multiworker_loader_preserves_seeded_shuffle_batches():
    dataset = fixed_dataset()
    reference = DataLoader(
        dataset,
        batch_size=8,
        shuffle=True,
        generator=seeded_generator(17),
        num_workers=0,
        pin_memory=False,
    )
    optimized = DataLoader(
        dataset,
        batch_size=8,
        shuffle=True,
        generator=seeded_generator(17),
        **optimized_loader_options(2),
    )

    for index, (left, right) in enumerate(zip(reference, optimized)):
        assert_nested_equal(left, right, f"batch[{index}]")


def test_nonblocking_cuda_transfer_preserves_prediction():
    if not torch.cuda.is_available():
        pytest.skip("CUDA equivalence requires the local drug environment")
    device = torch.device("cuda:0")
    model = torch.nn.Linear(8, 3, dtype=torch.float32).to(device).eval()
    features, labels = next(iter(DataLoader(fixed_dataset(), batch_size=16)))

    blocking = move_to_device((features, labels), device, non_blocking=False)
    pinned = features.pin_memory()
    nonblocking = move_to_device((pinned, labels.pin_memory()), device,
                                 non_blocking=True)
    torch.cuda.synchronize(device)

    assert_nested_equal(blocking, nonblocking)
    with torch.inference_mode():
        blocking_prediction = model(blocking[0])
        nonblocking_prediction = model(nonblocking[0])
    assert torch.equal(blocking_prediction, nonblocking_prediction)


def test_scalar_history_preserves_values_and_rejects_non_scalars():
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    values = [torch.tensor(value, dtype=torch.float64, device=device)
              for value in (1.25, -0.5, 3.0)]

    assert scalar_history_to_list(values) == [1.25, -0.5, 3.0]
    with pytest.raises(ValueError, match="scalar"):
        scalar_history_to_list([torch.ones(2, device=device)])


def test_deferred_history_preserves_training_trajectory_exactly():
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(23)
    reference = torch.nn.Linear(4, 1, dtype=torch.float64).to(device)
    deferred = deepcopy(reference)
    features = torch.arange(32, dtype=torch.float64, device=device).reshape(8, 4) / 32
    targets = torch.linspace(-0.25, 0.5, 8, dtype=torch.float64, device=device)[:, None]

    def train(model, defer):
        optimizer = torch.optim.SGD(model.parameters(), lr=0.05)
        history = []
        for _ in range(6):
            loss = (model(features) - targets).square().mean()
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()
            history.append(loss.detach() if defer else float(loss))
        if defer:
            history = scalar_history_to_list(history)
        return history, model(features).detach()

    reference_history, reference_prediction = train(reference, False)
    deferred_history, deferred_prediction = train(deferred, True)

    assert deferred_history == reference_history
    assert torch.equal(deferred_prediction, reference_prediction)
    for left, right in zip(reference.parameters(), deferred.parameters()):
        assert torch.equal(left, right)
