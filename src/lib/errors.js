/**
 * Structured error types for the BCDTPP library.
 */

/**
 * Thrown when mesh input data fails validation.
 */
export class MeshValidationError extends Error {
  /**
   * @param {string} message
   * @param {number=} index - The offending element index, if applicable.
   */
  constructor (message, index) {
    super(message)
    this.name = 'MeshValidationError'
    this.index = index
  }
}

/**
 * Thrown when a projection cannot be computed.
 */
export class ProjectionError extends Error {
  /**
   * @param {string} message - Description of the projection failure.
   */
  constructor (message) {
    super(message)
    this.name = 'ProjectionError'
  }
}

/**
 * Thrown when a linear system is singular or numerically ill-conditioned.
 */
export class SingularMatrixError extends Error {
  /**
   * @param {string} message - Description of the singularity condition.
   */
  constructor (message) {
    super(message)
    this.name = 'SingularMatrixError'
  }
}
