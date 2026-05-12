/**
 * Thrown when mesh input data fails validation.
 */
export class MeshValidationError extends Error {
  /**
   * @param message - Error message.
   * @param index - The offending element index, if applicable.
   */
  constructor(message: string, index?: number);
  name: 'MeshValidationError';
  /** The offending element index, if applicable. */
  index?: number;
}

/**
 * Thrown when a projection cannot be computed.
 */
export class ProjectionError extends Error {
  constructor(message: string);
  name: 'ProjectionError';
}

/**
 * Thrown when a linear system is singular or numerically ill-conditioned.
 */
export class SingularMatrixError extends Error {
  constructor(message: string);
  name: 'SingularMatrixError';
}
