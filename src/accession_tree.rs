use crate::accession_tree::AccessionTreeNode::{Accession, Branch};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub enum AccessionTreeNode {
    Accession(String),
    Branch(i32, i32),
}

#[derive(Serialize, Deserialize)]
pub struct AccessionTree {
    nodes: Vec<AccessionTreeNode>,
}

impl AccessionTree {
    pub fn new() -> Self {
        AccessionTree { nodes: Vec::new() }
    }

    /// Pushes the new node to the tree, returns the index of the new node
    pub fn push_new_node(&mut self, node: AccessionTreeNode) -> i32 {
        let index = self.nodes.len();
        self.nodes.push(node);
        index as i32
    }

    pub fn get_all_accession_indices(&self, index: i32) -> Vec<i32> {
        let mut vec = Vec::new();
        let mut current_node_index = index;
        while let Branch(left, right) = self.nodes.get(current_node_index as usize).unwrap() {
            vec.push(*right);
            current_node_index = *left;
        }
        vec.push(current_node_index);
        vec
    }

    pub fn get_accession_of_index(&self, accession_index: i32) -> &str {
        if let Accession(str) = self.nodes.get(accession_index as usize).unwrap() {
            &*str
        } else {
            panic!("Tried to look up an accession that wasn't an accession");
        }
    }
}
