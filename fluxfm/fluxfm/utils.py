class NotFittedError(ValueError, AttributeError):
    """Exception class to raise if estimator is used before fitting.
    References
    ----------
    .. [1] https://github.com/scikit-learn/scikit-learn/blob/0fb307bf39bbdacd6ed713c00724f8f871d60370/sklearn/exceptions.py#L19
    """

def check_is_fitted(estimator, attributes=None, *, msg=None, all_or_any=all):
    """Perform is_fitted validation for estimator.
    Checks if the estimator is fitted by verifying the presence of
    given fitted attributes (ending with a trailing underscore) and otherwise
    raises a NotFittedError with the given message.
    This utility is meant to be used internally by estimators themselves,
    typically in their own predict / transform methods.
    Parameters
    ----------
    estimator : estimator instance.
        estimator instance for which the check is performed.
    attributes : str, list or tuple of str, default=None
        Attribute name(s) given as string or a list/tuple of strings
        Eg.: ``["coef_", "estimator_", ...], "coef_"``
        If `None`, `estimator` is considered fitted if there exist an
        attribute that ends with a underscore and does not start with double
        underscore.
    msg : string
        The default error message is, "This %(name)s instance is not fitted
        yet. Call 'fit' with appropriate arguments before using this
        estimator."
        For custom messages if "%(name)s" is present in the message string,
        it is substituted for the estimator name.
        Eg. : "Estimator, %(name)s, must be fitted before sparsifying".
    all_or_any : callable, {all, any}, default all
        Specify whether all or any of the given attributes must exist.
    Returns
    -------
    None
    Raises
    ------
    NotFittedError
        If the attributes are not found.
    References
    ----------
    .. [1] https://github.com/scikit-learn/scikit-learn/blob/0fb307bf39bbdacd6ed713c00724f8f871d60370/sklearn/utils/validation.py#L956
    """
    if isclass(estimator):
        raise TypeError("{} is a class, not an instance.".format(estimator))
    if msg is None:
        msg = ("This %(name)s instance is not fitted yet. Call 'fit' with "
               "appropriate arguments before using this estimator.")

    if not hasattr(estimator, 'fit'):
        raise TypeError("%s is not an estimator instance." % (estimator))

    if attributes is not None:
        if not isinstance(attributes, (list, tuple)):
            attributes = [attributes]
        attrs = all_or_any([hasattr(estimator, attr) for attr in attributes])
    else:
        attrs = [v for v in vars(estimator)
                 if v.endswith("_") and not v.startswith("__")]

    if not attrs:
        raise NotFittedError(msg % {'name': type(estimator).__name__})
