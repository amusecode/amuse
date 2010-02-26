Datamodel
=========


All data is stored in sets
--------------------------
In the datamodel of AMUSE all data are stored in sets. The sets
store objects and the values of attributes belonging to the objects. 
All objects in a set have the same attributes, but not same values for
these attributes.

.. graphviz::
    
    digraph set0{
      fontsize=10.0;
      node [shape=record, fontsize=10.0,style=filled, fillcolor=lightyellow, width=3];
      struct1 [label="{set | object 1 | object 2 | .. | object n} | {attribute a | value 1.a| value 2.a | .. | value n.a} | {attribute b | value 1.b| value 2.b | .. | value n.b}| {.. | ..| .. | .. | ..} |  {attribute z | value 1.z| value 2.z | .. | value n.z}"];
    }
    

Like the relation database model
--------------------------------
For every object in a set, the set will store the values of the
attributes of the object. This model is like a relation 
database with Tables (sets in AMUSE), Records (an object in the set)
, Columns (an attribute of an object) and Fields (the value of an attribute of an object). 

.. graphviz::
    
    digraph layers0 {
      fontsize=10.0;
      node [fontsize=10.0,shape=box, style=filled, fillcolor=lightyellow, width=1.5];
      subgraph cluster0 {
      fontsize=10.0;
            style=filled;
            color=azure2;
            labeljust="l";
            label="AMUSE Model";
            "Set";
            "Object";
            "Attribute";
            "Value";
            "Set" -> "Attribute" [fontsize=10.0,label="defines"];
            "Object" -> "Value" [fontsize=10.0,label="has"];
            "Set" -> "Object" [fontsize=10.0,label="contains"];
            "Set" -> "Value"  [fontsize=10.0,label="stores"];
      } 
      subgraph cluster1 {
      fontsize=10.0;
            style=filled;
            color=azure2;
            labeljust="l";
            label="Relation Database Model";
            "Table";
            "Record";
            "Column";
            "Field";
            "Table" -> "Column" [fontsize=10.0,label="defines"];
            "Table" -> "Record" [fontsize=10.0,label="stores"];
            "Record" -> "Field" [fontsize=10.0,label="has"];
      }
    
    }

Objects are views
------------------
Objects from a set do not store any values, instead they defer
to the set to provide their attribute values. In a sense these
objects are pointers to a location in the set. When comparing
to the relational database model an object is like
a cursor. It can be used to access the values of the attributes
belonging to the object stored in the set.

.. graphviz::
    
    digraph layers0 {
      fontsize=10.0;
      node [fontsize=10.0, shape=box, style=filled, fillcolor=lightyellow, width=3];
      style=filled;
      color=azure2;
      "Set";
      "Object";
      "Object" -> "Set" [fontsize=10.0,label=" pointer to a location in"];
    }

When a user asks an object for its mass the object will query the
set for the stored value and return the answer to the user.

.. image:: objects_in_set.png

Objects have a key
------------------
The objects in a set can be identified with a unique key. All objects 
having the same key are seen as the same object by the AMUSE system.
The same object can exist in multiple sets. In each set this object
can have a different value of an attribute or different attributes. 
Or, in each set a different version of the object can exist. 

Sets use Storage Models
-----------------------
The actual storage of attribute values in a set is provided by a storage
model. The set provides the interface to the script writer, the 
storage model manages the data. Each storage model must decide how and
where to store the data. All data can be stored in the memory area
of the script or in the mememory area of the code or on a file or in
a relational database.

.. graphviz::
    
    digraph layers0 {
      fontsize=10.0;
      node [fontsize=10.0,shape=box, style=filled, fillcolor=lightyellow, width=1.5];
      subgraph cluster0 {
      fontsize=10.0;
            style=filled;
            color=azure2;
            labeljust="l";
            label="AMUSE Model";
            "Set";
      } 
      subgraph cluster1 {
      fontsize=10.0;
            style=filled;
            color=azure2;
            labeljust="l";
            label="Storage Models";
            "In Memory";
            "In Code";
            "In File";
            "In Database";
      }
      "Set" -> "In Memory";
      "Set" -> "In Code";
      "Set" -> "In File";
      "Set" -> "In Database";
    }


Selections on the set
---------------------
The datamodel provides subsets to handle a selection of the objects in a set.
When comparing to the relational database model an subset is like a view.
The subset does not store any data, all the data is stored in the original
set. When an attribute is updated in a subset, the attribute is also
updated in the original data.

 .. graphviz::
 
    digraph set0{
        fontsize=10.0;
        node [shape=record, fontsize=10.0,style=filled, fillcolor=lightyellow, width=4];
        subgraph cluster0 {
            fontsize=10.0;
            style=filled;
            color=azure2;
            labeljust="l";
            label="Subset";
            struct1 [label="{set |<here>  object 2 | .. | object m}"];
        }
 
        subgraph cluster1 {
            fontsize=10.0;
            style=filled;
            color=azure2;
            labeljust="l";
            label="Original Set";
            struct2 [label="{set | object 1 |<there> object 2 | .. | object n} | {attribute a | value 1.a| value 2.a | .. | value n.a} | {attribute b | value 1.b| value 2.b | .. | value n.b}| {.. | ..| .. | .. | ..} |  {attribute z | value 1.z| value 2.z | .. | value n.z}"];
        }
        
        struct1:here:e -> struct2:there:w
    
    }
    
