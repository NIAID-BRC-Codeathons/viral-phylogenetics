"""
Polyprotein Cleavage Site Prediction using ESM3 Embeddings

This module implements a binary classification approach for predicting viral
polyprotein cleavage sites. The approach has two main components:

1. ESM3 Embedding Generation: Extract per-residue embeddings from ESM3 for
   the entire polyprotein
2. Binary Classification: Feed embeddings to a small feedforward network to
   predict cut/no-cut

Author: AI Assistant
Date: November 4, 2025
"""

import torch
import torch.nn as nn
from torch.utils.data import Dataset
from typing import List, Tuple, Optional, Dict, Union
import numpy as np
import json
from dataclasses import dataclass
from tqdm import tqdm

# ESM imports
from esm.models.esm3 import ESM3
from esm.sdk.api import ESMProtein, LogitsConfig, ESM3InferenceClient


@dataclass
class CleavageDataPoint:
    """Data structure for a single polyprotein with cleavage annotations"""
    sequence: str
    cleavage_sites: List[int]  # 0-indexed positions where cleavage occurs
    protein_id: str
    organism: Optional[str] = None
    
    def get_binary_labels(self) -> List[int]:
        """Convert cleavage sites to binary labels (1=cut, 0=no-cut)"""
        labels = [0] * len(self.sequence)
        for site in self.cleavage_sites:
            if 0 <= site < len(self.sequence):
                labels[site] = 1
        return labels


class PolyproteinDataset(Dataset):
    """Dataset class for polyprotein cleavage prediction"""
    
    def __init__(self, data_points: List[CleavageDataPoint]):
        self.data_points = data_points
        
    def __len__(self):
        return len(self.data_points)
    
    def __getitem__(self, idx):
        data_point = self.data_points[idx]
        return {
            'sequence': data_point.sequence,
            'labels': torch.tensor(data_point.get_binary_labels(),
                                   dtype=torch.float32),
            'protein_id': data_point.protein_id
        }


class CleavagePredictionNetwork(nn.Module):
    """
    Small feedforward network for binary cleavage site prediction.
    Takes ESM3 embeddings as input and outputs binary predictions.
    """
    
    def __init__(self,
                 embedding_dim: int = 1536,  # ESM3-small embedding dimension
                 hidden_dims: List[int] = [512, 256, 128],
                 dropout_rate: float = 0.1):
        super().__init__()
        
        self.embedding_dim = embedding_dim
        
        # Build the feedforward layers
        layers = []
        prev_dim = embedding_dim
        
        for hidden_dim in hidden_dims:
            layers.extend([
                nn.Linear(prev_dim, hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout_rate)
            ])
            prev_dim = hidden_dim
        
        # Final classification layer
        layers.append(nn.Linear(prev_dim, 1))
        
        self.network = nn.Sequential(*layers)
        
        # Initialize weights
        self._initialize_weights()
    
    def _initialize_weights(self):
        """Initialize network weights using Xavier initialization"""
        for module in self.modules():
            if isinstance(module, nn.Linear):
                nn.init.xavier_uniform_(module.weight)
                nn.init.constant_(module.bias, 0.0)
    
    def forward(self, embeddings: torch.Tensor) -> torch.Tensor:
        """
        Forward pass
        
        Args:
            embeddings: Tensor of shape (batch_size, seq_len, embedding_dim)
            
        Returns:
            logits: Tensor of shape (batch_size, seq_len, 1)
        """
        batch_size, seq_len, embed_dim = embeddings.shape
        
        # Reshape to process all positions at once
        # (batch_size, seq_len, embed_dim) -> (batch_size * seq_len, embed_dim)
        embeddings_flat = embeddings.view(-1, embed_dim)
        
        # Forward through network
        logits_flat = self.network(embeddings_flat)
        
        # Reshape back to sequence format
        # (batch_size * seq_len, 1) -> (batch_size, seq_len, 1)
        logits = logits_flat.view(batch_size, seq_len, 1)
        
        return logits
    
    def predict_probabilities(self, embeddings: torch.Tensor) -> torch.Tensor:
        """Get cleavage probabilities using sigmoid activation"""
        logits = self.forward(embeddings)
        return torch.sigmoid(logits)


class ESM3EmbeddingExtractor:
    """Utility class for extracting embeddings from ESM3"""
    
    def __init__(self,
                 model: Union[ESM3, ESM3InferenceClient],
                 device: str = "cuda"):
        self.model = model
        self.device = device
        
        # Move model to device if it's a local model
        if isinstance(model, ESM3):
            self.model = self.model.to(device)
            self.model.eval()
    
    def extract_embeddings(self, sequence: str) -> torch.Tensor:
        """
        Extract per-residue embeddings for a protein sequence
        
        Args:
            sequence: Protein sequence string
            
        Returns:
            embeddings: Tensor of shape (seq_len, embedding_dim)
        """
        # Create ESMProtein object
        protein = ESMProtein(sequence=sequence)
        
        # Configure for embedding extraction
        config = LogitsConfig(return_embeddings=True)
        
        with torch.no_grad():
            if isinstance(self.model, ESM3):
                # Local model
                encoded = self.model.encode(protein)
                output = self.model.logits(encoded, config)
                embeddings = output.embeddings
            else:
                # Forge API client
                encoded = self.model.encode(protein)
                output = self.model.logits(encoded, config)
                embeddings = output.embeddings
        
        # Remove batch dimension and special tokens if needed
        if embeddings.dim() == 3:
            embeddings = embeddings.squeeze(0)  # Remove batch dim
        
        # Remove CLS and EOS tokens (first and last positions)
        if embeddings.size(0) == len(sequence) + 2:
            embeddings = embeddings[1:-1]  # Remove CLS and EOS
        elif embeddings.size(0) != len(sequence):
            raise ValueError(f"Embedding length {embeddings.size(0)} doesn't match sequence length {len(sequence)}")
        
        return embeddings.cpu()
    
    def extract_batch_embeddings(self, sequences: List[str]) -> List[torch.Tensor]:
        """Extract embeddings for a batch of sequences"""
        embeddings_list = []
        
        for sequence in tqdm(sequences, desc="Extracting embeddings"):
            embeddings = self.extract_embeddings(sequence)
            embeddings_list.append(embeddings)
        
        return embeddings_list


class PolyproteinCleavagePredictor:
    """
    Main class for training and inference of polyprotein cleavage site prediction
    """
    
    def __init__(self,
                 esm_model: Union[ESM3, ESM3InferenceClient],
                 device: str = "cuda",
                 embedding_dim: int = 1536):
        
        self.device = device
        self.embedding_extractor = ESM3EmbeddingExtractor(esm_model, device)
        self.classification_net = CleavagePredictionNetwork(embedding_dim=embedding_dim)
        self.classification_net.to(device)
        
        # Training history
        self.training_history = {
            'train_loss': [],
            'val_loss': [],
            'train_acc': [],
            'val_acc': []
        }
    
    def prepare_data(self, dataset: PolyproteinDataset) -> List[Dict]:
        """
        Prepare data by extracting embeddings and organizing with labels
        
        Args:
            dataset: PolyproteinDataset object
            
        Returns:
            List of dictionaries with embeddings and labels
        """
        prepared_data = []
        
        for i in tqdm(range(len(dataset)), desc="Preparing data"):
            item = dataset[i]
            sequence = item['sequence']
            labels = item['labels']
            protein_id = item['protein_id']
            
            # Extract embeddings
            embeddings = self.embedding_extractor.extract_embeddings(sequence)
            
            prepared_data.append({
                'embeddings': embeddings,
                'labels': labels,
                'protein_id': protein_id
            })
        
        return prepared_data
    
    def create_dataloader(self, prepared_data: List[Dict], batch_size: int = 8, shuffle: bool = True):
        """Create DataLoader for prepared data with proper collation"""
        
        def collate_fn(batch):
            """Custom collate function to handle variable-length sequences"""
            embeddings = [item['embeddings'] for item in batch]
            labels = [item['labels'] for item in batch]
            protein_ids = [item['protein_id'] for item in batch]
            
            # Pad sequences to same length
            max_len = max(emb.size(0) for emb in embeddings)
            embed_dim = embeddings[0].size(1)
            
            padded_embeddings = torch.zeros((len(batch), max_len, embed_dim))
            padded_labels = torch.zeros((len(batch), max_len))
            attention_mask = torch.zeros((len(batch), max_len))
            
            for i, (emb, lab) in enumerate(zip(embeddings, labels)):
                seq_len = emb.size(0)
                padded_embeddings[i, :seq_len] = emb
                padded_labels[i, :seq_len] = lab
                attention_mask[i, :seq_len] = 1.0
            
            return {
                'embeddings': padded_embeddings,
                'labels': padded_labels,
                'attention_mask': attention_mask,
                'protein_ids': protein_ids
            }
        
        from torch.utils.data import DataLoader
        
        return DataLoader(
            prepared_data,
            batch_size=batch_size,
            shuffle=shuffle,
            collate_fn=collate_fn,
            num_workers=0  # Set to 0 to avoid multiprocessing issues with ESM
        )
    
    def compute_metrics(self, predictions: torch.Tensor, labels: torch.Tensor, mask: torch.Tensor) -> Dict[str, float]:
        """Compute evaluation metrics"""
        # Apply mask to ignore padded positions
        valid_preds = predictions[mask.bool()]
        valid_labels = labels[mask.bool()]
        
        # Convert to binary predictions
        binary_preds = (valid_preds > 0.5).float()
        
        # Compute metrics
        correct = (binary_preds == valid_labels).sum().item()
        total = valid_labels.numel()
        accuracy = correct / total if total > 0 else 0.0
        
        # Compute precision, recall, F1 for positive class (cleavage sites)
        tp = ((binary_preds == 1) & (valid_labels == 1)).sum().item()
        fp = ((binary_preds == 1) & (valid_labels == 0)).sum().item()
        fn = ((binary_preds == 0) & (valid_labels == 1)).sum().item()
        
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0
        
        return {
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall,
            'f1_score': f1
        }
    
    def train_epoch(self, dataloader, optimizer, criterion):
        """Train for one epoch"""
        self.classification_net.train()
        total_loss = 0.0
        all_predictions = []
        all_labels = []
        all_masks = []
        
        for batch in tqdm(dataloader, desc="Training"):
            embeddings = batch['embeddings'].to(self.device)
            labels = batch['labels'].to(self.device)
            mask = batch['attention_mask'].to(self.device)
            
            optimizer.zero_grad()
            
            # Forward pass
            logits = self.classification_net(embeddings).squeeze(-1)  # Remove last dim
            predictions = torch.sigmoid(logits)
            
            # Compute loss only on valid positions
            loss = criterion(logits, labels)
            loss = (loss * mask).sum() / mask.sum()  # Masked average
            
            # Backward pass
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
            
            # Store for metrics computation
            all_predictions.append(predictions.detach().cpu())
            all_labels.append(labels.detach().cpu())
            all_masks.append(mask.detach().cpu())
        
        # Compute epoch metrics
        all_predictions = torch.cat(all_predictions, dim=0)
        all_labels = torch.cat(all_labels, dim=0)
        all_masks = torch.cat(all_masks, dim=0)
        
        metrics = self.compute_metrics(all_predictions, all_labels, all_masks)
        avg_loss = total_loss / len(dataloader)
        
        return avg_loss, metrics
    
    def validate_epoch(self, dataloader, criterion):
        """Validate for one epoch"""
        self.classification_net.eval()
        total_loss = 0.0
        all_predictions = []
        all_labels = []
        all_masks = []
        
        with torch.no_grad():
            for batch in tqdm(dataloader, desc="Validating"):
                embeddings = batch['embeddings'].to(self.device)
                labels = batch['labels'].to(self.device)
                mask = batch['attention_mask'].to(self.device)
                
                # Forward pass
                logits = self.classification_net(embeddings).squeeze(-1)
                predictions = torch.sigmoid(logits)
                
                # Compute loss
                loss = criterion(logits, labels)
                loss = (loss * mask).sum() / mask.sum()
                
                total_loss += loss.item()
                
                # Store for metrics computation
                all_predictions.append(predictions.detach().cpu())
                all_labels.append(labels.detach().cpu())
                all_masks.append(mask.detach().cpu())
        
        # Compute epoch metrics
        all_predictions = torch.cat(all_predictions, dim=0)
        all_labels = torch.cat(all_labels, dim=0)
        all_masks = torch.cat(all_masks, dim=0)
        
        metrics = self.compute_metrics(all_predictions, all_labels, all_masks)
        avg_loss = total_loss / len(dataloader)
        
        return avg_loss, metrics
    
    def train(self,
              train_dataset: PolyproteinDataset,
              val_dataset: Optional[PolyproteinDataset] = None,
              num_epochs: int = 50,
              batch_size: int = 8,
              learning_rate: float = 1e-3,
              weight_decay: float = 1e-5,
              save_path: Optional[str] = None):
        """
        Train the cleavage prediction model
        
        Args:
            train_dataset: Training dataset
            val_dataset: Validation dataset (optional)
            num_epochs: Number of training epochs
            batch_size: Batch size for training
            learning_rate: Learning rate for optimizer
            weight_decay: Weight decay for regularization
            save_path: Path to save the best model
        """
        
        # Prepare data
        print("Preparing training data...")
        train_data = self.prepare_data(train_dataset)
        train_loader = self.create_dataloader(train_data, batch_size, shuffle=True)
        
        if val_dataset is not None:
            print("Preparing validation data...")
            val_data = self.prepare_data(val_dataset)
            val_loader = self.create_dataloader(val_data, batch_size, shuffle=False)
        
        # Setup training
        optimizer = torch.optim.Adam(
            self.classification_net.parameters(),
            lr=learning_rate,
            weight_decay=weight_decay
        )
        
        # Use weighted BCE loss to handle class imbalance
        criterion = nn.BCEWithLogitsLoss(reduction='none')
        
        best_val_f1 = 0.0
        
        print(f"Starting training for {num_epochs} epochs...")
        
        for epoch in range(num_epochs):
            print(f"\nEpoch {epoch + 1}/{num_epochs}")
            
            # Training
            train_loss, train_metrics = self.train_epoch(train_loader, optimizer, criterion)
            self.training_history['train_loss'].append(train_loss)
            self.training_history['train_acc'].append(train_metrics['accuracy'])
            
            print(f"Train Loss: {train_loss:.4f}, Train Acc: {train_metrics['accuracy']:.4f}, "
                  f"Train F1: {train_metrics['f1_score']:.4f}")
            
            # Validation
            if val_dataset is not None:
                val_loss, val_metrics = self.validate_epoch(val_loader, criterion)
                self.training_history['val_loss'].append(val_loss)
                self.training_history['val_acc'].append(val_metrics['accuracy'])
                
                print(f"Val Loss: {val_loss:.4f}, Val Acc: {val_metrics['accuracy']:.4f}, "
                      f"Val F1: {val_metrics['f1_score']:.4f}")
                
                # Save best model
                if val_metrics['f1_score'] > best_val_f1:
                    best_val_f1 = val_metrics['f1_score']
                    if save_path:
                        self.save_model(save_path)
                        print(f"New best model saved with F1: {best_val_f1:.4f}")
    
    def predict(self, sequence: str) -> Tuple[np.ndarray, np.ndarray]:
        """
        Predict cleavage sites for a single sequence
        
        Args:
            sequence: Protein sequence string
            
        Returns:
            probabilities: Array of cleavage probabilities for each position
            binary_predictions: Array of binary predictions (0 or 1)
        """
        self.classification_net.eval()
        
        with torch.no_grad():
            # Extract embeddings
            embeddings = self.embedding_extractor.extract_embeddings(sequence)
            embeddings = embeddings.unsqueeze(0).to(self.device)  # Add batch dimension
            
            # Predict
            probabilities = self.classification_net.predict_probabilities(embeddings)
            probabilities = probabilities.squeeze().cpu().numpy()
            
            # Convert to binary predictions
            binary_predictions = (probabilities > 0.5).astype(int)
        
        return probabilities, binary_predictions
    
    def save_model(self, path: str):
        """Save the trained model"""
        torch.save({
            'model_state_dict': self.classification_net.state_dict(),
            'training_history': self.training_history,
        }, path)
        print(f"Model saved to {path}")
    
    def load_model(self, path: str):
        """Load a trained model"""
        checkpoint = torch.load(path, map_location=self.device)
        self.classification_net.load_state_dict(checkpoint['model_state_dict'])
        self.training_history = checkpoint.get('training_history', {})
        print(f"Model loaded from {path}")


def load_data_from_json(json_path: str) -> List[CleavageDataPoint]:
    """
    Load polyprotein data from JSON file
    
    Expected format:
    [
        {
            "protein_id": "protein_1",
            "sequence": "MAKL...",
            "cleavage_sites": [10, 25, 67],
            "organism": "SARS-CoV-2"
        },
        ...
    ]
    """
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    data_points = []
    for item in data:
        data_point = CleavageDataPoint(
            protein_id=item['protein_id'],
            sequence=item['sequence'],
            cleavage_sites=item['cleavage_sites'],
            organism=item.get('organism', None)
        )
        data_points.append(data_point)
    
    return data_points


# Example usage and training script
def main():
    """Example training script"""
    
    # Initialize ESM3 model
    print("Loading ESM3 model...")
    model = ESM3.from_pretrained("esm3_sm_open_v1").to("cuda")
    
    # Create predictor
    predictor = PolyproteinCleavagePredictor(model)
    
    # Example data (replace with your actual data loading)
    example_data = [
        CleavageDataPoint(
            protein_id="example_1",
            sequence="MKQHKAMIVALIVICITAVVAALVTRKDLCEVAKLR",
            cleavage_sites=[15, 25],  # Example cleavage positions
            organism="Example virus"
        )
    ]
    
    # Create datasets
    train_dataset = PolyproteinDataset(example_data)
    
    # Train model
    predictor.train(
        train_dataset=train_dataset,
        num_epochs=5,  # Use more epochs for real training
        batch_size=2,
        save_path="cleavage_model.pth"
    )
    
    # Make predictions
    test_sequence = "MKQHKAMIVALIVICITAVVAALVTRKDLCEVAKLR"
    probabilities, predictions = predictor.predict(test_sequence)
    
    print(f"Cleavage probabilities: {probabilities}")
    print(f"Binary predictions: {predictions}")


if __name__ == "__main__":
    main()