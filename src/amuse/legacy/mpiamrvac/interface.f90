MODULE mpiamrvac_interface

CONTAINS

    FUNCTION initialize_code()
        INTEGER initialize_code
        initialize_code = 0
    END FUNCTION
    
    FUNCTION cleanup_code()
        INTEGER cleanup_code
        cleanup_code = 0
    END FUNCTION
    
    FUNCTION commit_parameters()
        INTEGER commit_parameters
        commit_parameters = 0
    END FUNCTION
    
    FUNCTION recommit_parameters()
        INTEGER recommit_parameters
        recommit_parameters = 0
    END FUNCTION

END MODULE mpiamrvac_interface

